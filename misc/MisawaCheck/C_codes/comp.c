#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

int main(void){
 
  int            NPara_real,NPara_comp,NSlater;
  int            i;
  double         *Para_real,*Proj;
  double complex *Para_comp,*Slater,*OptTrans;

  NPara_comp = 20;
  Para_comp  = (double complex*)malloc(sizeof(double complex)*(NPara_comp));
  
  for(i=0;i<NPara_comp;i++){
    Para_comp[i] = 2*i+3*i*I;
  }
  
  for(i=0;i<NPara_comp;i++){
    printf("i=%d para = %lf %lf\n", i,&Para_comp[i],cimag(Para_comp[i])); 
  }

  //printf(" %d %d %d\n ", &Proj,&Slater,&OptTrans);

 return 0;
}
