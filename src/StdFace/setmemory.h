/*
HPhi-mVMC-StdFace - Common input generator
Copyright (C) 2015 The University of Tokyo

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
//
// Created by Kazuyoshi Yoshimi on 2019-01-09.
//

#ifndef MVMC_SETMEMORY_H
#define MVMC_SETMEMORY_H

#include <stdlib.h>
#include <complex.h>
///
/// \brief Allocation for A[N]
/// \param N [in] The size of the array A
/// \param A [in,out] Array to allocate
/// \return A
/// \author Kazuyoshi Yoshimi (University of Tokyo)
int *i_1d_allocate(const int N, int *A){
  A     = (int*)malloc((N)*sizeof(int));
  return A;
}

///
/// \brief Allocation for A[N][M]
/// \param N [in] The size of the array A
/// \param M [in] The size of the array M
/// \param A [in,out] Array to allocate
/// \return A
/// \author Kazuyoshi Yoshimi (University of Tokyo)
int **i_2d_allocate(int N, int M, int **A) {
  int int_i;
  A = (int **) malloc((N) * sizeof(int *));
  A[0] = (int *) malloc((M * N) * sizeof(int));
  for (int_i = 0; int_i < N; int_i++) {
    A[int_i] = A[0] + int_i * M;
  }
  return A;
}

///
/// \brief Allocation for A[N][M]
/// \param N [in] The size of the array A
/// \param M [in] The size of the array M
/// \param A [in,out] Array to allocate
/// \return A
/// \author Kazuyoshi Yoshimi (University of Tokyo)
double **d_2d_allocate(int N, int M, double **A){
  int int_i;
  A     = (double**)malloc((N)*sizeof(double*));
  A[0]  = (double*)malloc((M*N)*sizeof(double));
  for(int_i=0;int_i<N;int_i++){
    A[int_i] = A[0] + int_i*M;
  }
  return A;
}

///
/// \brief Allocation for A[N][M]
/// \param N [in] The size of the array A
/// \param M [in] The size of the array M
/// \param A [in,out] Array to allocate
/// \return A
/// \author Kazuyoshi Yoshimi (University of Tokyo)
complex double **cd_2d_allocate(int N, int M, complex double **A){
  int int_i;
  A     = (complex double**)malloc((N)*sizeof(complex double));
  A[0]  = (complex double*)malloc((M*N)*sizeof(complex double));
  for(int_i=0;int_i<N;int_i++){
    A[int_i] = A[0]+int_i*M;
  }
  return A;
}

//
/// \brief Allocation for A[N][M]
/// \param N [in] The size of the array A
/// \param M [in] The size of the array M
/// \param A [in,out] Array to allocate
/// \return A
/// \author Kazuyoshi Yoshimi (University of Tokyo)
double complex***cd_3d_allocate(int N, int M, int L, double complex***A){
  int int_i, int_j;
  A     = (double complex***)malloc((N)*sizeof(double complex**));
  A[0]  = (double complex**)malloc((M*N)*sizeof(double complex*));
  for(int_i=0;int_i<N; int_i++) {
    A[int_i] = A[0] + int_i*M;
    A[int_i][0] = (double complex*)malloc((M*L)*sizeof(double complex));
    for(int_j = 0; int_j<M; int_j++){
      A[int_i][int_j]= A[int_i][0] + int_j*L;
    }
  }
  return A;
}
#endif //MVMC_SETMEMORY_H
