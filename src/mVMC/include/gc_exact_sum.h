#ifndef MVMC_GC_EXACT_SUM_H
#define MVMC_GC_EXACT_SUM_H

#include <complex.h>

/* Exact zero of the sum of finite binary64 inputs, independently in real
 * and imaginary parts. Returns zero on success; invalid input returns -1
 * and clears *isZero. No floating-point additions or persistent storage. */
int GCExactComplexSumIsZero(const double complex *values, int count,
                          int *isZero);
/* Production initialization control: zero means the controls passed. */
int GCExactComplexSumSelfTest(void);

#endif
