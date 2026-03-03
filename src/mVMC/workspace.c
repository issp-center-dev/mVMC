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
/*-------------------------------------------------------------
 * Variational Monte Carlo
 * for workspace
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/
#include "workspace.h"
#pragma once

void initializeWorkSpaceAll() {
  int i;

  NumWorkSpaceInt = 0;
  WorkSpaceInt = NULL;

  NumWorkSpaceDouble = 0;
  WorkSpaceDouble = NULL;

  NumWorkSpaceComplex = 0; //TBC
  WorkSpaceComplex = NULL; //TBC

  NumWorkSpaceThreadInt = 0;
  WorkSpaceThreadInt = (int**)malloc(sizeof(int*)*NThread);
  WorkSpaceThreadIntNow = (int**)malloc(sizeof(int*)*NThread);

  NumWorkSpaceThreadDouble = 0;
  WorkSpaceThreadDouble = (double**)malloc(sizeof(double*)*NThread);
  WorkSpaceThreadDoubleNow = (double**)malloc(sizeof(double*)*NThread);

  NumWorkSpaceThreadComplex  = 0; // TBC
  WorkSpaceThreadComplex     = (double complex **)malloc(sizeof(double complex *)*NThread); //TBC
  WorkSpaceThreadComplexNow  = (double complex **)malloc(sizeof(double complex *)*NThread); //TBC

  for(i=0;i<NThread;i++) {
    WorkSpaceThreadInt[i] = NULL;
    WorkSpaceThreadDouble[i] = NULL;
    WorkSpaceThreadComplex[i] = NULL;
  }

  return;
}

void FreeWorkSpaceAll() {
  int i;

  for(i=0;i<NThread;i++) {
    free(WorkSpaceThreadInt[i]);
    free(WorkSpaceThreadDouble[i]);
    free(WorkSpaceThreadComplex[i]);
  }

  free(WorkSpaceThreadComplexNow);
  free(WorkSpaceThreadComplex);
  free(WorkSpaceThreadDoubleNow);
  free(WorkSpaceThreadDouble);
  free(WorkSpaceThreadIntNow);
  free(WorkSpaceThreadInt);
  free(WorkSpaceComplex);
  free(WorkSpaceDouble);
  free(WorkSpaceInt);

  return;
}

/* WorkSpaceInt */
void RequestWorkSpaceInt(int n) {
  if(n>NumWorkSpaceInt) {
    /* printf("WS Int %d -> %d\n",NumWorkSpaceInt,n); */
    NumWorkSpaceInt=n;
    free(WorkSpaceInt);
    WorkSpaceInt = (int*)malloc(sizeof(int)*NumWorkSpaceInt);
    WorkSpaceIntNow = WorkSpaceInt;
  }
  return;
}
int* GetWorkSpaceInt(int n) {
  int *p=WorkSpaceIntNow;
  WorkSpaceIntNow += n;
  return p;
}
void ReleaseWorkSpaceInt() {
  WorkSpaceIntNow = WorkSpaceInt;
  return;
}

/* WorkSpaceDouble */
void RequestWorkSpaceDouble(int n) {
  if(n>NumWorkSpaceDouble) {
    /* printf("WS Double %d -> %d\n",NumWorkSpaceDouble,n); */
    NumWorkSpaceDouble=n;
    free(WorkSpaceDouble);
    WorkSpaceDouble = (double*)malloc(sizeof(double)*NumWorkSpaceDouble);
    WorkSpaceDoubleNow = WorkSpaceDouble;
  }
  return;
}
double* GetWorkSpaceDouble(int n) {
  double *p=WorkSpaceDoubleNow;
  WorkSpaceDoubleNow += n;
  return p;
}
void ReleaseWorkSpaceDouble() {
  WorkSpaceDoubleNow = WorkSpaceDouble;
  return;
}

/* WorkSpaceComplex */
void RequestWorkSpaceComplex(int n) {
  if(n>NumWorkSpaceComplex) {
    NumWorkSpaceComplex=n;
    free(WorkSpaceComplex);
    WorkSpaceComplex = (double complex*)malloc(sizeof(double complex)*NumWorkSpaceComplex);
    WorkSpaceComplexNow = WorkSpaceComplex;
  }
  return;
}
double complex * GetWorkSpaceComplex(int n) {
  double complex *p=WorkSpaceComplexNow;
  WorkSpaceComplexNow += n;
  return p;
}
void ReleaseWorkSpaceComplex() {
  WorkSpaceComplexNow = WorkSpaceComplex;
  return;
}


/* WorkSpaceThreadInt */
void RequestWorkSpaceThreadInt(int n) {
  /* This function should be called "out of" an openMP parallel region. */
  int i;
  if(n>NumWorkSpaceThreadInt) {
    /* printf("WS ThreadInt %d -> %d\n",NumWorkSpaceThreadInt,n); */
    NumWorkSpaceThreadInt=n;
    #pragma omp parallel private(i)
    {
      i = omp_get_thread_num();
      free(WorkSpaceThreadInt[i]);
      WorkSpaceThreadInt[i] = (int*)malloc(sizeof(int)*NumWorkSpaceThreadInt);
      WorkSpaceThreadIntNow[i] = WorkSpaceThreadInt[i];
    }
  }
  return;
}
int* GetWorkSpaceThreadInt(int n) {
  /* This function should be called "in" an openMP parallel region. */
  int i = omp_get_thread_num();
  int *p = WorkSpaceThreadIntNow[i];
  WorkSpaceThreadIntNow[i] += n;
  return p;
}
void ReleaseWorkSpaceThreadInt() {
  /* This function should be called "out of" an openMP parallel region. */
  int i;
  for(i=0;i<NThread;i++) {
    WorkSpaceThreadIntNow[i] = WorkSpaceThreadInt[i];
  }
  return;
}

/* WorkSpaceThreadDouble */
void RequestWorkSpaceThreadDouble(int n) {
  /* This function should be called "out of" an openMP parallel region. */
  int i;
  if(n>NumWorkSpaceThreadDouble) {
    /* printf("WS ThreadDouble %d -> %d\n",NumWorkSpaceThreadDouble,n); */
    NumWorkSpaceThreadDouble=n;
    #pragma omp parallel private(i)
    {
      i = omp_get_thread_num();
      free(WorkSpaceThreadDouble[i]);
      WorkSpaceThreadDouble[i] = (double*)malloc(sizeof(double)*NumWorkSpaceThreadDouble);
      WorkSpaceThreadDoubleNow[i] = WorkSpaceThreadDouble[i];
    }
  }
  return;
}
double* GetWorkSpaceThreadDouble(int n) {
  /* This function should be called "in" an openMP parallel region. */
  int i = omp_get_thread_num();
  double *p = WorkSpaceThreadDoubleNow[i];
  WorkSpaceThreadDoubleNow[i] += n;
  return p;
}
void ReleaseWorkSpaceThreadDouble() {
  /* This function should be called "out of" an openMP parallel region. */
  int i;
  for(i=0;i<NThread;i++) {
    WorkSpaceThreadDoubleNow[i] = WorkSpaceThreadDouble[i];
  }
  return;
}

/* WorkSpaceThreadComplex */
void RequestWorkSpaceThreadComplex(int n) {
  /* This function should be called "out of" an openMP parallel region. */
  int i;
  if(n>NumWorkSpaceThreadComplex) {
    /* printf("WS ThreadComplex %d -> %d\n",NumWorkSpaceThreadComplex,n); */
    NumWorkSpaceThreadComplex=n;
    #pragma omp parallel private(i)
    {
      i = omp_get_thread_num();
      free(WorkSpaceThreadComplex[i]);
      WorkSpaceThreadComplex[i] = (double complex*)malloc(sizeof(double complex)*NumWorkSpaceThreadComplex);
      WorkSpaceThreadComplexNow[i] = WorkSpaceThreadComplex[i];
    }
  }
  return;
}
double complex* GetWorkSpaceThreadComplex(int n) {
  /* This function should be called "in" an openMP parallel region. */
  int i = omp_get_thread_num();
  double complex *p = WorkSpaceThreadComplexNow[i];
  WorkSpaceThreadComplexNow[i] += n;
  return p;
}
void ReleaseWorkSpaceThreadComplex() {
  /* This function should be called "out of" an openMP parallel region. */
  int i;
  for(i=0;i<NThread;i++) {
    WorkSpaceThreadComplexNow[i] = WorkSpaceThreadComplex[i];
  }
  return;
}
