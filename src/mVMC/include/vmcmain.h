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
 * main program header
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/

#ifndef _VMC_INCLUDE_FILES
#define _VMC_INCLUDE_FILES

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <errno.h>
#include <limits.h>
#include <string.h>
#include <complex.h>

#ifdef _mpi_use
  #include <mpi.h>
#else
typedef int MPI_Comm;
MPI_Comm MPI_COMM_WORLD=0;
inline void MPI_Init(int argc, char* argv[]) {return;}
inline void MPI_Finalize() {return;}
inline void MPI_Abort(MPI_Comm comm, int errorcode) {exit(errorcode); return;}
inline void MPI_Barrier(MPI_Comm comm) {return;}
inline void MPI_Comm_size(MPI_Comm comm, int *size) {*size = 1; return;}
inline void MPI_Comm_rank(MPI_Comm comm, int *rank) {*rank = 0; return;}
#endif /* _mpi_use */

extern int omp_get_max_threads(void);
extern int omp_get_thread_num(void);

#include "../../sfmt/SFMT.h"
#include "version.h"
#include "global.h"
#include "blas_externs.h"

#include "../safempi.c"
#include "../safempi_fcmp.c"
#include "../vmcclock.c"
#include "../workspace.c"

#ifdef _lapack
 #include "../stcopt_dposv.c"
#else
 #include "../stcopt_pdposv.c"
#endif

#include "../stcopt_cg.c"

#include "../stcopt.c"

#include "../gauleg.c"
#include "../legendrepoly.c"
#include "../avevar.c"
#include "../average.c"
#include "../parameter.c"
#include "../projection.c"
#include "../slater.c"
#include "../slater_fsz.c"
#include "../qp.c"
#include "../qp_real.c"
#include "../matrix.c"
#include "../pfupdate.c"
#include "../pfupdate_real.c"
#include "../pfupdate_fsz.c"
#include "../pfupdate_two_fcmp.c"
#include "../pfupdate_two_real.c"
#include "../pfupdate_two_fsz.c"
#include "../locgrn.c"
#include "../locgrn_real.c"
#include "../locgrn_fsz.c"
#include "../calham.c"
#include "../calham_real.c"
#include "../calham_fsz.c"
#include "../calgrn.c"
#include "../calgrn_fsz.c"
#include "../setmemory.c"
#include "../readdef.c"
#include "../initfile.c"

#include "../vmcmake.c"
#include "../vmcmake_real.c"
#include "../vmcmake_fsz.c"
#include "../vmccal.c"
#include "../vmccal_fsz.c"

#endif /* _VMC_INCLUDE_FILES */
