#ifndef _GC_ANTIPARALLEL_INPUT_H
#define _GC_ANTIPARALLEL_INPUT_H

#include <stdio.h>

/* Strict readers for the OrbitalAntiParallel (alias Orbital) file used by the
 * anti-parallel grand-canonical mode.  The caller owns fp; it is never closed
 * here.  Both return 0 on success and 1 on error, printing a diagnostic with
 * "GC anti-parallel", the file name and the line number to stderr. */

/* Read the five header lines from the current position: line 2 is
 * "label NOrbitalIdx" (> 0) and line 3 "label ComplexType" (0 or 1).  The
 * site count is checked so that (2*nsite)^2 fits in int.  Outputs are
 * written only on success. */
int GCAntiReadHeader(FILE *fp, int nsite, int *norb, int *complexFlag,
                     const char *filename);

/* Read the body right after the header: nsite*nsite rows "i j index [sign]"
 * (three or four integers; four with sign +1/-1 when apFlag is set, the
 * fourth column is otherwise ignored), then norb rows "index flag" in any
 * order, then only blank lines.  Every row is validated before its entries
 * are written.  OptFlag entries 2*(fidx+index) and 2*(fidx+index)+1 are set
 * and *optCount grows by one per row. */
int GCAntiReadOrbitals(FILE *fp, int **idx, int **sgn, int *opt,
                       int *optCount, int fidx, int complexFlag, int apFlag,
                       int nsite, int norb, const char *filename);

#endif
