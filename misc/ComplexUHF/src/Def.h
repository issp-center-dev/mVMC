#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <complex.h>
#define PI 3.14159265358979 
#define Fock 1
#include "struct.h"
#include "mfmemory.c"
#include "matrixlapack.c"
#include "readdef.c"
#include "check.c"
#include "initial.c"
#include "makeham.c"
#include "diag.c"
#include "green.c"
#include "cal_energy.c"
#include "output.c"
#include "cal_cisajs.c"
