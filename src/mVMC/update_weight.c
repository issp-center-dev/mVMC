/*
mVMC - A numerical solver package for a wide range of quantum lattice models based on many-variable Variational Monte Carlo method
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

#include "include/update_weight.h"

#define UPDATE_WEIGHT_LINE_MAX 512
#define UPDATE_WEIGHT_NAME_MAX 64

static const char *UpdateWeightNames[UpdateWeightCount] = {
    "Exchange", "LocalSpinFlip", "PairSpinFlip"};

static int SetUpdateWeightError(char *error, size_t errorSize,
                                const char *format, ...) {
  va_list arguments;
  if (error != NULL && errorSize > 0) {
    va_start(arguments, format);
    vsnprintf(error, errorSize, format, arguments);
    va_end(arguments);
  }
  return -1;
}

static int UpdateWeightNameIndex(const char *name) {
  int index;
  for (index = 0; index < UpdateWeightCount; index++) {
    if (strcmp(name, UpdateWeightNames[index]) == 0) return index;
  }
  return -1;
}

static int IsBlankOrComment(const char *line) {
  while (*line != '\0' && isspace((unsigned char)*line)) line++;
  return *line == '\0' || *line == '#';
}

int ReadUpdateWeight(FILE *fp, const char *source,
                     double weights[UpdateWeightCount],
                     char *error, size_t errorSize) {
  char line[UPDATE_WEIGHT_LINE_MAX];
  char keyword[UPDATE_WEIGHT_NAME_MAX];
  char extra[2];
  int count, index, row, scanned;
  int seen[UpdateWeightCount] = {0, 0, 0};
  double total = 0.0;

  if (fp == NULL || weights == NULL) {
    return SetUpdateWeightError(error, errorSize,
                                "invalid UpdateWeight parser arguments");
  }
  if (source == NULL) source = "UpdateWeight input";
  for (index = 0; index < UpdateWeightCount; index++) weights[index] = 0.0;

  if (fgets(line, sizeof(line), fp) == NULL ||
      fgets(line, sizeof(line), fp) == NULL) {
    return SetUpdateWeightError(error, errorSize,
                                "%s: incomplete UpdateWeight header", source);
  }
  scanned = sscanf(line, "%63s %d %1s", keyword, &count, extra);
  if (scanned != 2 || strcmp(keyword, "NUpdateWeight") != 0) {
    return SetUpdateWeightError(
        error, errorSize,
        "%s: expected 'NUpdateWeight <count>' on header line 2", source);
  }
  if (count < 1 || count > UpdateWeightCount) {
    return SetUpdateWeightError(
        error, errorSize, "%s: NUpdateWeight must be in [1, %d] (got %d)",
        source, UpdateWeightCount, count);
  }
  for (row = 0; row < 3; row++) {
    if (fgets(line, sizeof(line), fp) == NULL) {
      return SetUpdateWeightError(error, errorSize,
                                  "%s: incomplete five-line header", source);
    }
  }

  row = 0;
  while (fgets(line, sizeof(line), fp) != NULL) {
    double value;
    if (IsBlankOrComment(line)) continue;
    if (row >= count) {
      return SetUpdateWeightError(
          error, errorSize,
          "%s: more data rows than NUpdateWeight=%d", source, count);
    }
    scanned = sscanf(line, "%63s %lf %1s", keyword, &value, extra);
    if (scanned != 2) {
      return SetUpdateWeightError(
          error, errorSize, "%s: malformed UpdateWeight row %d", source,
          row + 1);
    }
    index = UpdateWeightNameIndex(keyword);
    if (index < 0) {
      return SetUpdateWeightError(
          error, errorSize, "%s: unknown update kernel '%s'", source,
          keyword);
    }
    if (seen[index] != 0) {
      return SetUpdateWeightError(
          error, errorSize, "%s: duplicate update kernel '%s'", source,
          keyword);
    }
    if (!isfinite(value) || value < 0.0) {
      return SetUpdateWeightError(
          error, errorSize,
          "%s: weight for %s must be finite and non-negative", source,
          keyword);
    }
    seen[index] = 1;
    weights[index] = value;
    total += value;
    row++;
  }
  if (row != count) {
    return SetUpdateWeightError(
        error, errorSize, "%s: expected %d data rows but read %d", source,
        count, row);
  }
  if (!isfinite(total) || total <= 0.0) {
    return SetUpdateWeightError(
        error, errorSize, "%s: total update weight must be positive", source);
  }
  for (index = 0; index < UpdateWeightCount; index++) {
    weights[index] /= total;
  }
  return 0;
}

int ValidateUpdateWeightContract(
    int enabled, int updatePath, int twoSz, int nsite, int nLocSpn, int ne,
    int orbitalGeneral, const double weights[UpdateWeightCount],
    char *error, size_t errorSize) {
  if (!enabled) return 0;
  if (weights == NULL) {
    return SetUpdateWeightError(error, errorSize,
                                "UpdateWeight values are unavailable");
  }
  if (updatePath != 2) {
    return SetUpdateWeightError(
        error, errorSize,
        "InUpdateWeight currently requires NExUpdatePath=2 (got %d)",
        updatePath);
  }
  if (weights[UpdateWeightLocalSpinFlip] > 0.0 &&
      (twoSz != -1 || orbitalGeneral == 0)) {
    return SetUpdateWeightError(
        error, errorSize,
        "LocalSpinFlip with NExUpdatePath=2 requires 2Sz=-1 and "
        "OrbitalGeneral");
  }
  if (weights[UpdateWeightPairSpinFlip] <= 0.0) return 0;
  if (twoSz != -1) {
    return SetUpdateWeightError(
        error, errorSize, "PairSpinFlip requires 2Sz=-1 (got %d)", twoSz);
  }
  if (orbitalGeneral == 0) {
    return SetUpdateWeightError(error, errorSize,
                                "PairSpinFlip requires OrbitalGeneral");
  }
  if (nsite < 2 || nLocSpn != nsite || 2LL * ne != nsite) {
    return SetUpdateWeightError(
        error, errorSize,
        "PairSpinFlip requires Nsite>=2, NLocalSpin=Nsite, and Ncond=0 "
        "(Nsite=%d, NLocalSpin=%d, 2*Ne=%lld)",
        nsite, nLocSpn, 2LL * ne);
  }
  return 0;
}

int ValidateUpdateWeightLocSpin(
    int enabled, int nsite, int nLocSpn, const int *locSpn,
    const double weights[UpdateWeightCount], char *error, size_t errorSize) {
  int site;
  if (!enabled || weights == NULL ||
      weights[UpdateWeightPairSpinFlip] <= 0.0) return 0;
  if (locSpn == NULL || nLocSpn != nsite) {
    return SetUpdateWeightError(
        error, errorSize,
        "PairSpinFlip requires one local-spin entry for every site");
  }
  for (site = 0; site < nsite; site++) {
    if (locSpn[site] != 1) {
      return SetUpdateWeightError(
          error, errorSize,
          "PairSpinFlip requires LocSpin=1 at every site (site %d has %d)",
          site, locSpn[site]);
    }
  }
  return 0;
}

int SelectUpdateWeight(const double weights[UpdateWeightCount], double draw) {
  double cumulative = 0.0;
  int index, lastPositive = -1;
  if (weights == NULL || !isfinite(draw) || draw < 0.0 || draw >= 1.0)
    return -1;
  for (index = 0; index < UpdateWeightCount; index++) {
    if (weights[index] > 0.0) lastPositive = index;
    cumulative += weights[index];
    if (draw < cumulative) return index;
  }
  return lastPositive;
}
