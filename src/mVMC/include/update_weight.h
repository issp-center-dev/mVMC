/*
mVMC - A numerical solver package for a wide range of quantum lattice models based on many-variable Variational Monte Carlo method
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/
#pragma once

#include <stddef.h>
#include <stdio.h>

enum UpdateWeightIndex {
  UpdateWeightExchange = 0,
  UpdateWeightLocalSpinFlip,
  UpdateWeightPairSpinFlip,
  UpdateWeightCount
};

int ReadUpdateWeight(FILE *fp, const char *source,
                     double weights[UpdateWeightCount],
                     char *error, size_t errorSize);
int ValidateUpdateWeightContract(
    int enabled, int updatePath, int twoSz, int nsite, int nLocSpn, int ne,
    int orbitalGeneral, const double weights[UpdateWeightCount],
    char *error, size_t errorSize);
int ValidateUpdateWeightLocSpin(
    int enabled, int nsite, int nLocSpn, const int *locSpn,
    const double weights[UpdateWeightCount], char *error, size_t errorSize);
int SelectUpdateWeight(const double weights[UpdateWeightCount], double draw);
