#ifndef _VMCCAL_FSZ
#define _VMCCAL_FSZ
#include <complex.h>

/* parent: run-wide diagnostics/accounting; child: sample decomposition. */
void VMCMainCal_fsz(MPI_Comm comm_parent, MPI_Comm comm_child);
void VMC_BF_MainCal_fsz(MPI_Comm comm_parent, MPI_Comm comm_child);

#endif
