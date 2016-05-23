#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

int main(void){
 
  int            NPara_real,NPara_comp,NSlater;
  int            i;
  double         *Para_real,*Proj;
  double complex *Para_comp,*Slater,*OptTrans;

  NPara_real = 10;
  NPara_comp = 23;
  NSlater    = 5;

  Para_real = (double*)malloc(sizeof(double)*(NPara_real));
  Para_comp = (double complex*)malloc(sizeof(double complex)*(NPara_comp));
  
  for(i=0;i<NPara_real;i++){
    Para_real[i] = i;
  }
  for(i=0;i<NPara_comp;i++){
    Para_comp[i] = 2*i+3*i*I;
  }
  
  Proj     = Para_real;
  Slater   = Para_comp;
  OptTrans = Para_comp + NSlater;

  for(i=0;i<NPara_real+1;i++){
    printf("i=%d para = %lf \n", i,Proj[i]); 
  }
  printf(" \n"); 
  for(i=0;i<NSlater;i++){
    printf("i=%d para = %lf %lf\n", i,creal(Slater[i]),cimag(Slater[i])); 
  }
  printf(" \n"); 
  for(i=0;i<NPara_comp-NSlater;i++){
    printf("i=%d para = %lf %lf\n", i,creal(OptTrans[i]),cimag(OptTrans[i])); 
  }

  //printf(" %d %d %d\n ", &Proj,&Slater,&OptTrans);

 return 0;
}
