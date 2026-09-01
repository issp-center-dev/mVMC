#ifndef _MATRIX_GC_H
#define _MATRIX_GC_H

#define GC_MALL_OK 0
#define GC_MALL_PFAFFIAN 1
#define GC_MALL_GETRF 2
#define GC_MALL_GETRI 3
#define GC_MALL_NONFINITE 4
#define GC_MALL_INVALID_ARGUMENT 5

int CalculateMAllGC_fcmp(int ncur, const int *eleIdx, int qpStart,
                         int qpEnd);

#endif
