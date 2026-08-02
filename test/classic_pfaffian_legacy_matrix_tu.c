#include <complex.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#ifdef _mpi_use
#include <mpi.h>
#else
typedef int MPI_Comm;
#endif

extern int omp_get_thread_num(void);

#include "global.h"
#include "blas_externs.h"

/* Compile the legacy unity-only matrix implementation as an isolated oracle. */
#include "../src/mVMC/matrix.c"

/* Compile the actual classic rank-one/two fast proposal kernels as observers. */
#include "../src/mVMC/pfupdate.c"
#include "../src/mVMC/pfupdate_real.c"
#include "../src/mVMC/pfupdate_two_fcmp.c"
#include "../src/mVMC/pfupdate_two_real.c"
