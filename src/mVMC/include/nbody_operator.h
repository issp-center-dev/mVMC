#ifndef MVMC_NBODY_OPERATOR_H
#define MVMC_NBODY_OPERATOR_H

typedef enum {
  NBODY_REDUCED_ZERO = 0,
  NBODY_REDUCED_SCALAR = 1,
  NBODY_REDUCED_HOPS = 2,
  NBODY_REDUCED_INVALID = 3
} NBodyReducedKind;

typedef struct {
  NBodyReducedKind kind;
  int order;
  int sign;
} NBodyReduction;

int ReduceNBodyTerm(int n, const int *rsi, const int *rsj,
                    const int *occupation, int nOrbitals,
                    int *rsiWork, int *rsjWork, int workCapacity,
                    NBodyReduction *reduction);

#endif
