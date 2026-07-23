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

/*
 * Reduce an N-body operator without modifying the caller's input arrays.
 *
 * rsi and rsj each contain n orbital indices, occupation contains nOrbitals
 * binary entries, and workCapacity is the capacity (in int elements) of both
 * rsiWork and rsjWork.  The writable ranges rsiWork[0:n], rsjWork[0:n], and
 * *reduction must not overlap each other or any readable input range.
 * Read-only input ranges may overlap each other.
 *
 * On a valid call, reduction is fully initialized.  A physical-zero result
 * may leave the work arrays partially reduced; callers must use reduction.kind
 * rather than inspect unused work entries.  Invalid non-aliasing calls set
 * reduction to NBODY_REDUCED_INVALID before returning -1.  Calls that alias
 * reduction with an array are rejected without writing through reduction.
 */
int ReduceNBodyTerm(int n, const int *rsi, const int *rsj,
                    const int *occupation, int nOrbitals,
                    int *rsiWork, int *rsjWork, int workCapacity,
                    NBodyReduction *reduction);

#endif
