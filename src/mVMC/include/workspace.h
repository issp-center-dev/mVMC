#ifndef _WORKSPACE
#define _WORKSPACE
#include <complex.h>
int NumWorkSpaceInt;
int *WorkSpaceInt;
int *WorkSpaceIntNow;

int NumWorkSpaceDouble;
double *WorkSpaceDouble;
double *WorkSpaceDoubleNow;

int NumWorkSpaceComplex;
double complex *WorkSpaceComplex;
double complex *WorkSpaceComplexNow;

int NumWorkSpaceThreadInt;
int **WorkSpaceThreadInt;
int **WorkSpaceThreadIntNow;

int NumWorkSpaceThreadDouble;
double **WorkSpaceThreadDouble;
double **WorkSpaceThreadDoubleNow;

int NumWorkSpaceThreadComplex;
double complex **WorkSpaceThreadComplex;
double complex **WorkSpaceThreadComplexNow;

void initializeWorkSpaceAll();
void FreeWorkSpaceAll();

void RequestWorkSpaceInt(int n);
int* GetWorkSpaceInt(int n);
void ReleaseWorkSpaceInt();

void RequestWorkSpaceDouble(int n);
double* GetWorkSpaceDouble(int n);
void ReleaseWorkSpaceDouble();

void RequestWorkSpaceComplex(int n);
double complex* GetWorkSpaceComplex(int n);
void ReleaseWorkSpaceComplex();

void RequestWorkSpaceThreadInt(int n);
int* GetWorkSpaceThreadInt(int n);
void ReleaseWorkSpaceThreadInt();

void RequestWorkSpaceThreadDouble(int n);
double* GetWorkSpaceThreadDouble(int n);
void ReleaseWorkSpaceThreadDouble();

void requestWorkSpaceThreadComplex(int n);
double complex* GetWorkSpaceThreadComplex(int n);
void ReleaseWorkSpaceThreadComplex();

#endif
