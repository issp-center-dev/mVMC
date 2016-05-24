#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

int main(void){
 
  int            NPara_real;
  int            i,j;
  int         *tmp,**Para_real,*pDouble;

  NPara_real = 2;
  tmp        = (int*)malloc(sizeof(int)*(NPara_real*NPara_real));
  //pDouble    = tmp; 
  Para_real  = (int**)malloc(sizeof(int*)*(NPara_real));

  for(i=0;i<NPara_real;i++){
    Para_real[i]  = tmp+i*NPara_real;
    printf("i=%d: tmp=%08X  %08X %08X\n",i,tmp,tmp+i*NPara_real,Para_real[i]);
    //pDouble      += NPara_real;
  }
  for(i=0;i<NPara_real;i++){
    for(j=0;j<NPara_real;j++){
      printf("%d %d %08X\n",i,j,&(Para_real[i][j]));
      Para_real[i][j]= i*j;
    } 
  }
  for(i=0;i<NPara_real;i++){
    for(j=0;j<NPara_real;j++){
      printf("%i=d j=%d %d \n",i,j,*Para_real[i]);
      printf("%d %d %d \n",i,j,Para_real[i][j]);
    } 
  }

 return 0;
}
