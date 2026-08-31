/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#include "classic_krylov_model.h"

#if defined(MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE) ||                     \
    defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

struct MVMCClassicKrylovModelWorkspace {
  MVMCKrylovFockModel model;
  MVMCKrylovHamiltonianTerm *terms;
  MVMCKrylovFermionOperator *operators;
  size_t allocated_bytes;
};

typedef struct {
  int values[12];
} RawMetadata;

typedef struct {
  MVMCClassicKrylovTransfer *transfers;
  MVMCClassicKrylovSiteCoupling *coulomb_intra;
  MVMCClassicKrylovPairCoupling *coulomb_inter;
  MVMCClassicKrylovPairCoupling *hund;
  MVMCClassicKrylovPairCoupling *exchange;
} RawCopies;

static int checked_multiply(size_t left, size_t right, size_t *product) {
  if (product == NULL || (left != 0 && right > SIZE_MAX / left)) return 0;
  *product = left * right;
  return 1;
}

static int checked_add(size_t left, size_t right, size_t *sum) {
  if (sum == NULL || right > SIZE_MAX - left) return 0;
  *sum = left + right;
  return 1;
}

static int valid_site(int site, int site_count) {
  return site >= 0 && site < site_count;
}

static int valid_real(double value) {
  return isfinite(value);
}

static int valid_complex(double complex value) {
  return isfinite(creal(value)) && isfinite(cimag(value));
}

static MVMCKrylovStatus validate_root_raw(
    const MVMCClassicKrylovRawModel *raw) {
  int spin_changing_transfer = 0;
  int index;

  if (raw == NULL || raw->site_count <= 0 ||
      raw->site_count > INT_MAX / 2 || raw->up_electron_count <= 0 ||
      raw->up_electron_count != raw->down_electron_count ||
      raw->up_electron_count > raw->site_count ||
      (raw->pure_spin != 0 && raw->pure_spin != 1) ||
      (raw->pure_spin &&
       raw->up_electron_count + raw->down_electron_count !=
           raw->site_count) ||
      raw->transfer_count < 0 || raw->coulomb_intra_count < 0 ||
      raw->coulomb_inter_count < 0 || raw->hund_count < 0 ||
      raw->exchange_count < 0 || raw->pair_hopping_count < 0 ||
      raw->inter_all_count < 0 || raw->nbody_inter_all_count < 0) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if ((raw->transfer_count != 0 && raw->transfers == NULL) ||
      (raw->coulomb_intra_count != 0 && raw->coulomb_intra == NULL) ||
      (raw->coulomb_inter_count != 0 && raw->coulomb_inter == NULL) ||
      (raw->hund_count != 0 && raw->hund == NULL) ||
      (raw->exchange_count != 0 && raw->exchange == NULL)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  for (index = 0; index < raw->transfer_count; ++index) {
    const MVMCClassicKrylovTransfer *transfer = raw->transfers + index;
    if (!valid_site(transfer->output_site, raw->site_count) ||
        !valid_site(transfer->input_site, raw->site_count) ||
        (transfer->output_spin != 0 && transfer->output_spin != 1) ||
        (transfer->input_spin != 0 && transfer->input_spin != 1) ||
        !valid_complex(transfer->coefficient)) {
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    }
    if (transfer->output_spin != transfer->input_spin) {
      spin_changing_transfer = 1;
    }
  }
  for (index = 0; index < raw->coulomb_intra_count; ++index) {
    if (!valid_site(raw->coulomb_intra[index].site, raw->site_count) ||
        !valid_real(raw->coulomb_intra[index].coefficient)) {
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    }
  }
#define VALIDATE_PAIRS(field, count)                                      \
  do {                                                                     \
    for (index = 0; index < (count); ++index) {                            \
      if (!valid_site((field)[index].first_site, raw->site_count) ||       \
          !valid_site((field)[index].second_site, raw->site_count) ||      \
          !valid_real((field)[index].coefficient)) {                       \
        return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;                        \
      }                                                                    \
    }                                                                      \
  } while (0)
  VALIDATE_PAIRS(raw->coulomb_inter, raw->coulomb_inter_count);
  VALIDATE_PAIRS(raw->hund, raw->hund_count);
  VALIDATE_PAIRS(raw->exchange, raw->exchange_count);
#undef VALIDATE_PAIRS
  if (spin_changing_transfer || raw->pair_hopping_count != 0 ||
      raw->inter_all_count != 0 || raw->nbody_inter_all_count != 0 ||
      (!raw->pure_spin && raw->exchange_count != 0) ||
      (raw->pure_spin &&
       (raw->transfer_count != 0 || raw->coulomb_intra_count != 0))) {
    return MVMC_KRYLOV_STATUS_UNSUPPORTED_MODEL;
  }
  return MVMC_KRYLOV_STATUS_OK;
}

static void metadata_from_raw(const MVMCClassicKrylovRawModel *raw,
                              RawMetadata *metadata) {
  metadata->values[0] = raw->site_count;
  metadata->values[1] = raw->up_electron_count;
  metadata->values[2] = raw->down_electron_count;
  metadata->values[3] = raw->pure_spin;
  metadata->values[4] = raw->transfer_count;
  metadata->values[5] = raw->coulomb_intra_count;
  metadata->values[6] = raw->coulomb_inter_count;
  metadata->values[7] = raw->hund_count;
  metadata->values[8] = raw->exchange_count;
  metadata->values[9] = raw->pair_hopping_count;
  metadata->values[10] = raw->inter_all_count;
  metadata->values[11] = raw->nbody_inter_all_count;
}

static void raw_from_metadata(const RawMetadata *metadata,
                              MVMCClassicKrylovRawModel *raw,
                              const RawCopies *copies) {
  memset(raw, 0, sizeof(*raw));
  raw->site_count = metadata->values[0];
  raw->up_electron_count = metadata->values[1];
  raw->down_electron_count = metadata->values[2];
  raw->pure_spin = metadata->values[3];
  raw->transfer_count = metadata->values[4];
  raw->transfers = copies->transfers;
  raw->coulomb_intra_count = metadata->values[5];
  raw->coulomb_intra = copies->coulomb_intra;
  raw->coulomb_inter_count = metadata->values[6];
  raw->coulomb_inter = copies->coulomb_inter;
  raw->hund_count = metadata->values[7];
  raw->hund = copies->hund;
  raw->exchange_count = metadata->values[8];
  raw->exchange = copies->exchange;
  raw->pair_hopping_count = metadata->values[9];
  raw->inter_all_count = metadata->values[10];
  raw->nbody_inter_all_count = metadata->values[11];
}

static int allocate_array(int count, size_t item_size, void **array) {
  size_t bytes;
  *array = NULL;
  if (count == 0) return 1;
  if (count < 0 || !checked_multiply((size_t)count, item_size, &bytes)) {
    return 0;
  }
  *array = malloc(bytes);
  return *array != NULL;
}

static int allocate_raw_copies(const RawMetadata *metadata,
                               RawCopies *copies) {
  memset(copies, 0, sizeof(*copies));
  return allocate_array(metadata->values[4], sizeof(*copies->transfers),
                        (void **)&copies->transfers) &&
         allocate_array(metadata->values[5], sizeof(*copies->coulomb_intra),
                        (void **)&copies->coulomb_intra) &&
         allocate_array(metadata->values[6], sizeof(*copies->coulomb_inter),
                        (void **)&copies->coulomb_inter) &&
         allocate_array(metadata->values[7], sizeof(*copies->hund),
                        (void **)&copies->hund) &&
         allocate_array(metadata->values[8], sizeof(*copies->exchange),
                        (void **)&copies->exchange);
}

static void free_raw_copies(RawCopies *copies) {
  if (copies == NULL) return;
  free(copies->exchange);
  free(copies->hund);
  free(copies->coulomb_inter);
  free(copies->coulomb_intra);
  free(copies->transfers);
  memset(copies, 0, sizeof(*copies));
}

#ifdef _mpi_use
static int broadcast_bytes(void *data, size_t bytes,
                           MVMCClassicPfaffianCommunicator communicator) {
  unsigned char *position = (unsigned char *)data;
  size_t offset = 0;

  while (offset < bytes) {
    const size_t remaining = bytes - offset;
    const int count = remaining > (size_t)INT_MAX ? INT_MAX : (int)remaining;
    if (MPI_Bcast(position + offset, count, MPI_BYTE, 0, communicator) !=
        MPI_SUCCESS) {
      return 0;
    }
    offset += (size_t)count;
  }
  return 1;
}
#endif

static int copy_or_broadcast_array(void *destination, const void *root_source,
                                   int count, size_t item_size,
                                   int rank,
                                   MVMCClassicPfaffianCommunicator communicator) {
  size_t bytes;
  if (count == 0) return 1;
  if (!checked_multiply((size_t)count, item_size, &bytes)) return 0;
  if (destination == NULL) return 0;
  if (rank == 0 && root_source == NULL) return 0;
  if (rank == 0) memcpy(destination, root_source, bytes);
#ifdef _mpi_use
  return broadcast_bytes(destination, bytes, communicator);
#else
  (void)communicator;
  return 1;
#endif
}

static int transfer_key_compare(const void *left_pointer,
                                const void *right_pointer) {
  const MVMCClassicKrylovTransfer *left =
      (const MVMCClassicKrylovTransfer *)left_pointer;
  const MVMCClassicKrylovTransfer *right =
      (const MVMCClassicKrylovTransfer *)right_pointer;
#define CMP_FIELD(field)                  \
  do {                                    \
    if (left->field < right->field) return -1; \
    if (left->field > right->field) return 1;  \
  } while (0)
  CMP_FIELD(output_site);
  CMP_FIELD(output_spin);
  CMP_FIELD(input_site);
  CMP_FIELD(input_spin);
#undef CMP_FIELD
  if (creal(left->coefficient) < creal(right->coefficient)) return -1;
  if (creal(left->coefficient) > creal(right->coefficient)) return 1;
  if (cimag(left->coefficient) < cimag(right->coefficient)) return -1;
  if (cimag(left->coefficient) > cimag(right->coefficient)) return 1;
  return 0;
}

static int transfer_same_key(const MVMCClassicKrylovTransfer *left,
                             const MVMCClassicKrylovTransfer *right) {
  return left->output_site == right->output_site &&
         left->output_spin == right->output_spin &&
         left->input_site == right->input_site &&
         left->input_spin == right->input_spin;
}

static int site_compare(const void *left_pointer, const void *right_pointer) {
  const MVMCClassicKrylovSiteCoupling *left =
      (const MVMCClassicKrylovSiteCoupling *)left_pointer;
  const MVMCClassicKrylovSiteCoupling *right =
      (const MVMCClassicKrylovSiteCoupling *)right_pointer;
  if (left->site < right->site) return -1;
  if (left->site > right->site) return 1;
  if (left->coefficient < right->coefficient) return -1;
  if (left->coefficient > right->coefficient) return 1;
  return 0;
}

static int pair_compare(const void *left_pointer, const void *right_pointer) {
  const MVMCClassicKrylovPairCoupling *left =
      (const MVMCClassicKrylovPairCoupling *)left_pointer;
  const MVMCClassicKrylovPairCoupling *right =
      (const MVMCClassicKrylovPairCoupling *)right_pointer;
  if (left->first_site < right->first_site) return -1;
  if (left->first_site > right->first_site) return 1;
  if (left->second_site < right->second_site) return -1;
  if (left->second_site > right->second_site) return 1;
  if (left->coefficient < right->coefficient) return -1;
  if (left->coefficient > right->coefficient) return 1;
  return 0;
}

static int add_real(double value, double *sum, double *compensation) {
  const double updated = *sum + value;
  if (!isfinite(value) || !isfinite(updated)) return 0;
  if (fabs(*sum) >= fabs(value)) {
    *compensation += (*sum - updated) + value;
  } else {
    *compensation += (value - updated) + *sum;
  }
  if (!isfinite(*compensation)) return 0;
  *sum = updated;
  return 1;
}

static MVMCKrylovStatus compact_transfers(MVMCClassicKrylovTransfer *values,
                                          int *count) {
  int input = 0;
  int output = 0;
  if (*count == 0) return MVMC_KRYLOV_STATUS_OK;
  qsort(values, (size_t)*count, sizeof(*values), transfer_key_compare);
  while (input < *count) {
    const MVMCClassicKrylovTransfer key = values[input];
    double real_sum = 0.0, real_compensation = 0.0;
    double imag_sum = 0.0, imag_compensation = 0.0;
    do {
      if (!add_real(creal(values[input].coefficient), &real_sum,
                    &real_compensation) ||
          !add_real(cimag(values[input].coefficient), &imag_sum,
                    &imag_compensation)) {
        return MVMC_KRYLOV_STATUS_NONFINITE;
      }
      ++input;
    } while (input < *count && transfer_same_key(&key, values + input));
    real_sum += real_compensation;
    imag_sum += imag_compensation;
    if (!isfinite(real_sum) || !isfinite(imag_sum)) {
      return MVMC_KRYLOV_STATUS_NONFINITE;
    }
    if (real_sum != 0.0 || imag_sum != 0.0) {
      values[output] = key;
      values[output].coefficient = real_sum + I * imag_sum;
      ++output;
    }
  }
  *count = output;
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus compact_sites(MVMCClassicKrylovSiteCoupling *values,
                                      int *count) {
  int input = 0;
  int output = 0;
  if (*count == 0) return MVMC_KRYLOV_STATUS_OK;
  qsort(values, (size_t)*count, sizeof(*values), site_compare);
  while (input < *count) {
    const int site = values[input].site;
    double sum = 0.0, compensation = 0.0;
    do {
      if (!add_real(values[input].coefficient, &sum, &compensation)) {
        return MVMC_KRYLOV_STATUS_NONFINITE;
      }
      ++input;
    } while (input < *count && values[input].site == site);
    sum += compensation;
    if (!isfinite(sum)) return MVMC_KRYLOV_STATUS_NONFINITE;
    if (sum != 0.0) {
      values[output].site = site;
      values[output].coefficient = sum;
      ++output;
    }
  }
  *count = output;
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus compact_pairs(MVMCClassicKrylovPairCoupling *values,
                                      int *count) {
  int index;
  int input = 0;
  int output = 0;
  if (*count == 0) return MVMC_KRYLOV_STATUS_OK;
  for (index = 0; index < *count; ++index) {
    if (values[index].second_site < values[index].first_site) {
      const int temporary = values[index].first_site;
      values[index].first_site = values[index].second_site;
      values[index].second_site = temporary;
    }
  }
  qsort(values, (size_t)*count, sizeof(*values), pair_compare);
  while (input < *count) {
    const int first = values[input].first_site;
    const int second = values[input].second_site;
    double sum = 0.0, compensation = 0.0;
    do {
      if (!add_real(values[input].coefficient, &sum, &compensation)) {
        return MVMC_KRYLOV_STATUS_NONFINITE;
      }
      ++input;
    } while (input < *count && values[input].first_site == first &&
             values[input].second_site == second);
    sum += compensation;
    if (!isfinite(sum)) return MVMC_KRYLOV_STATUS_NONFINITE;
    if (sum != 0.0) {
      values[output].first_site = first;
      values[output].second_site = second;
      values[output].coefficient = sum;
      ++output;
    }
  }
  *count = output;
  return MVMC_KRYLOV_STATUS_OK;
}

static const MVMCClassicKrylovTransfer *find_transfer(
    const MVMCClassicKrylovTransfer *values, int count,
    int output_site, int output_spin, int input_site, int input_spin) {
  int lower = 0;
  int upper = count;
  while (lower < upper) {
    const int middle = lower + (upper - lower) / 2;
    const MVMCClassicKrylovTransfer *candidate = values + middle;
    int relation = 0;
    if (candidate->output_site != output_site) {
      relation = candidate->output_site < output_site ? -1 : 1;
    } else if (candidate->output_spin != output_spin) {
      relation = candidate->output_spin < output_spin ? -1 : 1;
    } else if (candidate->input_site != input_site) {
      relation = candidate->input_site < input_site ? -1 : 1;
    } else if (candidate->input_spin != input_spin) {
      relation = candidate->input_spin < input_spin ? -1 : 1;
    } else {
      return candidate;
    }
    if (relation < 0) {
      lower = middle + 1;
    } else {
      upper = middle;
    }
  }
  return NULL;
}

static MVMCKrylovStatus validate_hermitian_transfers(
    const MVMCClassicKrylovTransfer *values, int count) {
  int index;
  for (index = 0; index < count; ++index) {
    const MVMCClassicKrylovTransfer *value = values + index;
    if (value->output_site == value->input_site &&
        value->output_spin == value->input_spin) {
      if (cimag(value->coefficient) != 0.0) {
        return MVMC_KRYLOV_STATUS_NON_HERMITIAN_MODEL;
      }
    } else {
      const MVMCClassicKrylovTransfer *reverse = find_transfer(
          values, count, value->input_site, value->input_spin,
          value->output_site, value->output_spin);
      if (reverse == NULL) return MVMC_KRYLOV_STATUS_NON_HERMITIAN_MODEL;
      if (value->coefficient != conj(reverse->coefficient)) {
        return MVMC_KRYLOV_STATUS_NON_HERMITIAN_MODEL;
      }
    }
  }
  return MVMC_KRYLOV_STATUS_OK;
}

static int append_term(MVMCClassicKrylovModelWorkspace *workspace,
                       size_t *term_index, size_t *operator_index,
                       size_t term_capacity, size_t operator_capacity,
                       double complex coefficient,
                       MVMCClassicKrylovSourceKind source_kind,
                       size_t source_index,
                       const size_t *orbitals, size_t orbital_count) {
  MVMCKrylovHamiltonianTerm *term;
  size_t index;
  if (orbital_count != 2 && orbital_count != 4) return 0;
  if (workspace == NULL || workspace->terms == NULL ||
      workspace->operators == NULL || term_index == NULL ||
      operator_index == NULL || *term_index >= term_capacity ||
      *operator_index > operator_capacity ||
      orbital_count > operator_capacity - *operator_index) {
    return 0;
  }
  term = workspace->terms + *term_index;
  term->coefficient = coefficient;
  term->operator_offset = *operator_index;
  term->operator_count = orbital_count;
  term->source_kind = (int)source_kind;
  term->source_index = source_index;
  for (index = 0; index < orbital_count; ++index) {
    workspace->operators[*operator_index + index].kind =
        (index % 2 == 0) ? MVMC_KRYLOV_FERMION_CREATE
                         : MVMC_KRYLOV_FERMION_ANNIHILATE;
    workspace->operators[*operator_index + index].orbital = orbitals[index];
  }
  ++*term_index;
  *operator_index += orbital_count;
  return 1;
}

static int append_density(MVMCClassicKrylovModelWorkspace *workspace,
                          size_t *term_index, size_t *operator_index,
                          size_t term_capacity, size_t operator_capacity,
                          double coefficient,
                          MVMCClassicKrylovSourceKind source_kind,
                          size_t source_index,
                          size_t first_orbital, size_t second_orbital) {
  const size_t orbitals[4] = {first_orbital, first_orbital,
                              second_orbital, second_orbital};
  return append_term(workspace, term_index, operator_index,
                     term_capacity, operator_capacity, coefficient,
                     source_kind, source_index, orbitals, 4);
}

static MVMCKrylovStatus build_model(MVMCClassicKrylovRawModel *raw,
                                    MVMCClassicKrylovModelWorkspace **output) {
  MVMCClassicKrylovModelWorkspace *workspace = NULL;
  MVMCKrylovStatus status;
  size_t term_count = 0;
  size_t operator_count = 0;
  size_t term_bytes = 0;
  size_t operator_bytes = 0;
  size_t term_index = 0;
  size_t operator_index = 0;
  int transfer_count = raw->transfer_count;
  int intra_count = raw->coulomb_intra_count;
  int inter_count = raw->coulomb_inter_count;
  int hund_count = raw->hund_count;
  int exchange_count = raw->exchange_count;
  int index;

  *output = NULL;
  status = compact_transfers((MVMCClassicKrylovTransfer *)raw->transfers,
                             &transfer_count);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = compact_sites((MVMCClassicKrylovSiteCoupling *)raw->coulomb_intra,
                           &intra_count);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = compact_pairs((MVMCClassicKrylovPairCoupling *)raw->coulomb_inter,
                           &inter_count);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = compact_pairs((MVMCClassicKrylovPairCoupling *)raw->hund,
                           &hund_count);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = compact_pairs((MVMCClassicKrylovPairCoupling *)raw->exchange,
                           &exchange_count);
  }
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  status = validate_hermitian_transfers(raw->transfers, transfer_count);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;

#define ADD_TERMS(count, multiplier)                                      \
  do {                                                                     \
    size_t addition;                                                       \
    if (!checked_multiply((size_t)(count), (size_t)(multiplier),           \
                          &addition) ||                                    \
        !checked_add(term_count, addition, &term_count)) {                 \
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;                          \
    }                                                                      \
  } while (0)
  ADD_TERMS(transfer_count, 1);
  ADD_TERMS(intra_count, 1);
  ADD_TERMS(inter_count, 4);
  ADD_TERMS(hund_count, 2);
  ADD_TERMS(exchange_count, 2);
#undef ADD_TERMS
  if (!checked_multiply((size_t)transfer_count, 2, &operator_count)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
#define ADD_OPERATORS(count, multiplier)                                  \
  do {                                                                     \
    size_t addition;                                                       \
    if (!checked_multiply((size_t)(count), (size_t)(multiplier),           \
                          &addition) ||                                    \
        !checked_add(operator_count, addition, &operator_count)) {         \
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;                          \
    }                                                                      \
  } while (0)
  ADD_OPERATORS(intra_count, 4);
  ADD_OPERATORS(inter_count, 16);
  ADD_OPERATORS(hund_count, 8);
  ADD_OPERATORS(exchange_count, 8);
#undef ADD_OPERATORS
  if (!checked_multiply(term_count, sizeof(*workspace->terms), &term_bytes) ||
      !checked_multiply(operator_count, sizeof(*workspace->operators),
                        &operator_bytes)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  workspace = (MVMCClassicKrylovModelWorkspace *)calloc(1, sizeof(*workspace));
  if (workspace != NULL && term_count != 0) {
    workspace->terms =
        (MVMCKrylovHamiltonianTerm *)calloc(term_count,
                                            sizeof(*workspace->terms));
  }
  if (workspace != NULL && operator_count != 0) {
    workspace->operators =
        (MVMCKrylovFermionOperator *)calloc(operator_count,
                                             sizeof(*workspace->operators));
  }
  if (workspace == NULL || (term_count != 0 && workspace->terms == NULL) ||
      (operator_count != 0 && workspace->operators == NULL)) {
    mvmc_classic_krylov_model_workspace_destroy(workspace);
    return MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
  }
#define APPEND_OR_FAIL(expression)                                        \
  do {                                                                    \
    if (!(expression)) {                                                  \
      mvmc_classic_krylov_model_workspace_destroy(workspace);             \
      return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;               \
    }                                                                     \
  } while (0)
  for (index = 0; index < transfer_count; ++index) {
    const MVMCClassicKrylovTransfer *transfer = raw->transfers + index;
    const size_t orbitals[2] = {
        (size_t)transfer->output_site +
            (size_t)transfer->output_spin * (size_t)raw->site_count,
        (size_t)transfer->input_site +
            (size_t)transfer->input_spin * (size_t)raw->site_count};
    APPEND_OR_FAIL(append_term(
        workspace, &term_index, &operator_index, term_count, operator_count,
        -transfer->coefficient, MVMC_CLASSIC_KRYLOV_SOURCE_TRANSFER,
        (size_t)index, orbitals, 2));
  }
  for (index = 0; index < intra_count; ++index) {
    const size_t site = (size_t)raw->coulomb_intra[index].site;
    APPEND_OR_FAIL(append_density(
        workspace, &term_index, &operator_index, term_count, operator_count,
        raw->coulomb_intra[index].coefficient,
        MVMC_CLASSIC_KRYLOV_SOURCE_COULOMB_INTRA, (size_t)index, site,
        site + (size_t)raw->site_count));
  }
  for (index = 0; index < inter_count; ++index) {
    int first_spin, second_spin;
    for (first_spin = 0; first_spin < 2; ++first_spin) {
      for (second_spin = 0; second_spin < 2; ++second_spin) {
        APPEND_OR_FAIL(append_density(
            workspace, &term_index, &operator_index,
            term_count, operator_count,
            raw->coulomb_inter[index].coefficient,
            MVMC_CLASSIC_KRYLOV_SOURCE_COULOMB_INTER,
            (size_t)index * 4 + (size_t)(2 * first_spin + second_spin),
            (size_t)raw->coulomb_inter[index].first_site +
                (size_t)first_spin * (size_t)raw->site_count,
            (size_t)raw->coulomb_inter[index].second_site +
                (size_t)second_spin * (size_t)raw->site_count));
      }
    }
  }
  for (index = 0; index < hund_count; ++index) {
    int spin;
    for (spin = 0; spin < 2; ++spin) {
      APPEND_OR_FAIL(append_density(
          workspace, &term_index, &operator_index,
          term_count, operator_count,
          -raw->hund[index].coefficient,
          MVMC_CLASSIC_KRYLOV_SOURCE_HUND,
          (size_t)index * 2 + (size_t)spin,
          (size_t)raw->hund[index].first_site +
              (size_t)spin * (size_t)raw->site_count,
          (size_t)raw->hund[index].second_site +
              (size_t)spin * (size_t)raw->site_count));
    }
  }
  for (index = 0; index < exchange_count; ++index) {
    const size_t first = (size_t)raw->exchange[index].first_site;
    const size_t second = (size_t)raw->exchange[index].second_site;
    const size_t up_first = first;
    const size_t up_second = second;
    const size_t down_first = first + (size_t)raw->site_count;
    const size_t down_second = second + (size_t)raw->site_count;
    const size_t forward[4] = {up_first, up_second,
                               down_second, down_first};
    const size_t reverse[4] = {up_second, up_first,
                               down_first, down_second};
    APPEND_OR_FAIL(append_term(
        workspace, &term_index, &operator_index, term_count, operator_count,
        raw->exchange[index].coefficient,
        MVMC_CLASSIC_KRYLOV_SOURCE_EXCHANGE,
        (size_t)index * 2, forward, 4));
    APPEND_OR_FAIL(append_term(
        workspace, &term_index, &operator_index, term_count, operator_count,
        raw->exchange[index].coefficient,
        MVMC_CLASSIC_KRYLOV_SOURCE_EXCHANGE,
        (size_t)index * 2 + 1, reverse, 4));
  }
#undef APPEND_OR_FAIL
  if (term_index != term_count || operator_index != operator_count) {
    mvmc_classic_krylov_model_workspace_destroy(workspace);
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  workspace->model.site_count = (size_t)raw->site_count;
  workspace->model.up_electron_count = (size_t)raw->up_electron_count;
  workspace->model.down_electron_count = (size_t)raw->down_electron_count;
  workspace->model.pure_spin = raw->pure_spin;
  workspace->model.hermitian = 1;
  workspace->model.terms = workspace->terms;
  workspace->model.term_count = term_count;
  workspace->model.operators = workspace->operators;
  workspace->model.operator_count = operator_count;
  if (!checked_add(sizeof(*workspace), term_bytes,
                   &workspace->allocated_bytes) ||
      !checked_add(workspace->allocated_bytes, operator_bytes,
                   &workspace->allocated_bytes)) {
    mvmc_classic_krylov_model_workspace_destroy(workspace);
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  *output = workspace;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_classic_krylov_model_workspace_create_from_root(
    const MVMCClassicKrylovRawModel *raw,
    MVMCClassicPfaffianCommunicator communicator,
    MVMCClassicKrylovModelWorkspace **workspace) {
  MVMCClassicKrylovRawModel local_raw;
  MVMCClassicKrylovModelWorkspace *created = NULL;
  RawMetadata metadata;
  RawCopies copies;
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_OK;
  int rank = 0;
  int local_ready;
  int all_ready;
  int local_status;
  int global_status;
  int broadcast_ok = 1;

  if (workspace == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  *workspace = NULL;
  memset(&metadata, 0, sizeof(metadata));
  memset(&copies, 0, sizeof(copies));
#ifdef _mpi_use
  if (communicator == MPI_COMM_NULL ||
      MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
#else
  (void)communicator;
#endif
  if (rank == 0) {
    status = validate_root_raw(raw);
    if (status == MVMC_KRYLOV_STATUS_OK) metadata_from_raw(raw, &metadata);
  }
  local_status = (int)status;
#ifdef _mpi_use
  if (MPI_Bcast(&local_status, 1, MPI_INT, 0, communicator) != MPI_SUCCESS ||
      MPI_Bcast(&metadata, (int)sizeof(metadata), MPI_BYTE, 0, communicator) !=
          MPI_SUCCESS) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
#endif
  status = (MVMCKrylovStatus)local_status;
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  local_ready = allocate_raw_copies(&metadata, &copies);
#ifdef _mpi_use
  if (MPI_Allreduce(&local_ready, &all_ready, 1, MPI_INT, MPI_MIN,
                    communicator) != MPI_SUCCESS) {
    free_raw_copies(&copies);
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
#else
  all_ready = local_ready;
#endif
  if (!all_ready) {
    free_raw_copies(&copies);
    return MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
  }
#define COPY_FIELD(field, count_index)                                    \
  do {                                                                     \
    const void *root_source = rank == 0 ? (const void *)raw->field : NULL; \
    if (!copy_or_broadcast_array(                                          \
            copies.field, root_source, metadata.values[count_index],       \
            sizeof(*copies.field), rank, communicator)) {                  \
      broadcast_ok = 0;                                                    \
    }                                                                      \
  } while (0)
  COPY_FIELD(transfers, 4);
  COPY_FIELD(coulomb_intra, 5);
  COPY_FIELD(coulomb_inter, 6);
  COPY_FIELD(hund, 7);
  COPY_FIELD(exchange, 8);
#undef COPY_FIELD
#ifdef _mpi_use
  if (MPI_Allreduce(&broadcast_ok, &all_ready, 1, MPI_INT, MPI_MIN,
                    communicator) != MPI_SUCCESS) {
    free_raw_copies(&copies);
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
#else
  all_ready = broadcast_ok;
#endif
  if (!all_ready) {
    free_raw_copies(&copies);
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  raw_from_metadata(&metadata, &local_raw, &copies);
  status = build_model(&local_raw, &created);
  local_status = (int)status;
#ifdef _mpi_use
  if (MPI_Allreduce(&local_status, &global_status, 1, MPI_INT, MPI_MAX,
                    communicator) != MPI_SUCCESS) {
    mvmc_classic_krylov_model_workspace_destroy(created);
    free_raw_copies(&copies);
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
#else
  global_status = local_status;
#endif
  free_raw_copies(&copies);
  if (global_status != (int)MVMC_KRYLOV_STATUS_OK) {
    mvmc_classic_krylov_model_workspace_destroy(created);
    return (MVMCKrylovStatus)global_status;
  }
  *workspace = created;
  return MVMC_KRYLOV_STATUS_OK;
}

void mvmc_classic_krylov_model_workspace_destroy(
    MVMCClassicKrylovModelWorkspace *workspace) {
  if (workspace == NULL) return;
  free(workspace->operators);
  free(workspace->terms);
  free(workspace);
}

const MVMCKrylovFockModel *mvmc_classic_krylov_model(
    const MVMCClassicKrylovModelWorkspace *workspace) {
  return workspace == NULL ? NULL : &workspace->model;
}

size_t mvmc_classic_krylov_model_workspace_bytes(
    const MVMCClassicKrylovModelWorkspace *workspace) {
  return workspace == NULL ? 0 : workspace->allocated_bytes;
}

#endif /* reference or bounded engine */
