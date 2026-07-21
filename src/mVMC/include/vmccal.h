#ifndef _VMCCAL
#define _VMCCAL
#include <complex.h>
/* parent: run-wide diagnostics/accounting; child: sample decomposition. */
void VMCMainCal(MPI_Comm comm_parent, MPI_Comm comm_child);
void VMC_BF_MainCal(MPI_Comm comm_parent, MPI_Comm comm_child);
#endif
