#ifndef _BACKFLOW_STABLE
#define _BACKFLOW_STABLE

#include <complex.h>

/* Scratch for the non-FSZ NSPGaussLeg==1 occupied kernel.  The caller owns
 * every buffer so proposal and Green-function loops do not allocate memory. */
typedef struct {
  double *factor;               /* Ne*Ne */
  double *work;                 /* lwork; required when inverse != NULL */
  int *iwork;                   /* Nsize; row map, then LU pivots */
  int lwork;
} BFStableWorkspaceReal;

typedef struct {
  double complex *factor;       /* Ne*Ne */
  double complex *work;         /* lwork; required when inverse != NULL */
  int *iwork;                   /* Nsize; row map, then LU pivots */
  int lwork;
} BFStableWorkspaceFcmp;

#endif
