/*-------------------------------------------------------------
 * Variational Monte Carlo
 * for workspace
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/

int NumWorkSpaceInt;
int *WorkSpaceInt;
int *WorkSpaceIntNow;

int NumWorkSpaceDouble;
double *WorkSpaceDouble;
double *WorkSpaceDoubleNow;

int NumWorkSpaceThreadInt;
int **WorkSpaceThreadInt;
int **WorkSpaceThreadIntNow;

int NumWorkSpaceThreadDouble;
double **WorkSpaceThreadDouble;
double **WorkSpaceThreadDoubleNow;

void initializeWorkSpaceAll();
void FreeWorkSpaceAll();

void RequestWorkSpaceInt(int n);
int* GetWorkSpaceInt(int n);
void ReleaseWorkSpaceInt();

void RequestWorkSpaceDouble(int n);
double* GetWorkSpaceDouble(int n);
void ReleaseWorkSpaceDouble();

void RequestWorkSpaceThreadInt(int n);
int* GetWorkSpaceThreadInt(int n);
void ReleaseWorkSpaceThreadInt();

void RequestWorkSpaceThreadDouble(int n);
double* GetWorkSpaceThreadDouble(int n);
void ReleaseWorkSpaceThreadDouble();

void initializeWorkSpaceAll() {
  int i;

  NumWorkSpaceInt = 0;
  WorkSpaceInt = NULL;

  NumWorkSpaceDouble = 0;
  WorkSpaceDouble = NULL;

  NumWorkSpaceThreadInt = 0;
  WorkSpaceThreadInt = (int**)malloc(sizeof(int*)*NThread);
  WorkSpaceThreadIntNow = (int**)malloc(sizeof(int*)*NThread);

  NumWorkSpaceThreadDouble = 0;
  WorkSpaceThreadDouble = (double**)malloc(sizeof(double*)*NThread);
  WorkSpaceThreadDoubleNow = (double**)malloc(sizeof(double*)*NThread);

  for(i=0;i<NThread;i++) {
    WorkSpaceThreadInt[i] = NULL;
    WorkSpaceThreadDouble[i] = NULL;
  }

  return;
}

void FreeWorkSpaceAll() {
  int i;

  for(i=0;i<NThread;i++) {
    free(WorkSpaceThreadInt[i]);
    free(WorkSpaceThreadDouble[i]);
  }

  free(WorkSpaceThreadDoubleNow);
  free(WorkSpaceThreadDouble);
  free(WorkSpaceThreadIntNow);
  free(WorkSpaceThreadInt);
  free(WorkSpaceDouble);
  free(WorkSpaceInt);

  return;
}

/* WorkSpaceInt */
void RequestWorkSpaceInt(int n) {
  if(n>NumWorkSpaceInt) {
    /* printf("WS Int %d -> %d\n",NumWorkSpaceInt,n); */
    NumWorkSpaceInt=n;
    free(WorkSpaceInt);
    WorkSpaceInt = (int*)malloc(sizeof(int)*NumWorkSpaceInt);
    WorkSpaceIntNow = WorkSpaceInt;
  }
  return;
}
int* GetWorkSpaceInt(int n) {
  int *p=WorkSpaceIntNow;
  WorkSpaceIntNow += n;
  return p;
}
void ReleaseWorkSpaceInt() {
  WorkSpaceIntNow = WorkSpaceInt;
  return;
}

/* WorkSpaceDouble */
void RequestWorkSpaceDouble(int n) {
  if(n>NumWorkSpaceDouble) {
    /* printf("WS Double %d -> %d\n",NumWorkSpaceDouble,n); */
    NumWorkSpaceDouble=n;
    free(WorkSpaceDouble);
    WorkSpaceDouble = (double*)malloc(sizeof(double)*NumWorkSpaceDouble);
    WorkSpaceDoubleNow = WorkSpaceDouble;
  }
  return;
}
double* GetWorkSpaceDouble(int n) {
  double *p=WorkSpaceDoubleNow;
  WorkSpaceDoubleNow += n;
  return p;
}
void ReleaseWorkSpaceDouble() {
  WorkSpaceDoubleNow = WorkSpaceDouble;
  return;
}

/* WorkSpaceThreadInt */
void RequestWorkSpaceThreadInt(int n) {
  /* This function should be called "out of" an openMP parallel region. */
  int i;
  if(n>NumWorkSpaceThreadInt) {
    /* printf("WS ThreadInt %d -> %d\n",NumWorkSpaceThreadInt,n); */
    NumWorkSpaceThreadInt=n;
    #pragma omp parallel private(i)
    {
      i = omp_get_thread_num();
      free(WorkSpaceThreadInt[i]);
      WorkSpaceThreadInt[i] = (int*)malloc(sizeof(int)*NumWorkSpaceThreadInt);
      WorkSpaceThreadIntNow[i] = WorkSpaceThreadInt[i];
    }
  }
  return;
}
int* GetWorkSpaceThreadInt(int n) {
  /* This function should be called "in" an openMP parallel region. */
  int i = omp_get_thread_num();
  int *p = WorkSpaceThreadIntNow[i];
  WorkSpaceThreadIntNow[i] += n;
  return p;
}
void ReleaseWorkSpaceThreadInt() {
  /* This function should be called "out of" an openMP parallel region. */
  int i;
  for(i=0;i<NThread;i++) {
    WorkSpaceThreadIntNow[i] = WorkSpaceThreadInt[i];
  }
  return;
}

/* WorkSpaceThreadDouble */
void RequestWorkSpaceThreadDouble(int n) {
  /* This function should be called "out of" an openMP parallel region. */
  int i;
  if(n>NumWorkSpaceThreadDouble) {
    /* printf("WS ThreadDouble %d -> %d\n",NumWorkSpaceThreadDouble,n); */
    NumWorkSpaceThreadDouble=n;
    #pragma omp parallel private(i)
    {
      i = omp_get_thread_num();
      free(WorkSpaceThreadDouble[i]);
      WorkSpaceThreadDouble[i] = (double*)malloc(sizeof(double)*NumWorkSpaceThreadDouble);
      WorkSpaceThreadDoubleNow[i] = WorkSpaceThreadDouble[i];
    }
  }
  return;
}
double* GetWorkSpaceThreadDouble(int n) {
  /* This function should be called "in" an openMP parallel region. */
  int i = omp_get_thread_num();
  double *p = WorkSpaceThreadDoubleNow[i];
  WorkSpaceThreadDoubleNow[i] += n;
  return p;
}
void ReleaseWorkSpaceThreadDouble() {
  /* This function should be called "out of" an openMP parallel region. */
  int i;
  for(i=0;i<NThread;i++) {
    WorkSpaceThreadDoubleNow[i] = WorkSpaceThreadDouble[i];
  }
  return;
}
