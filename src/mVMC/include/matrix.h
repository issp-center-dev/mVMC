#ifndef _MATRIX
#define _MATRIX

int CalculateMAll_fcmp(const int *eleIdx, const int qpStart, const int qpEnd);

int CalculateMAll_real(const int *eleIdx, const int qpStart, const int qpEnd);

int CalculateMAll_BF_real(const int *eleIdx, const int qpStart, const int qpEnd);

int CalculateMAll_BF_fcmp(const int *eleIdx, const int qpStart, const int qpEnd);

#endif


