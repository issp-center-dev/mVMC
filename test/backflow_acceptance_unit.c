#include <float.h>
#include <stdio.h>
#include "../src/mVMC/include/backflow_acceptance.h"

int main(void) {
  const struct { double p, n, o, u; int valid, accept; } cases[]={
    {0,500,0,0.999,1,1}, /* finite log ratio +1000 must accept */
    {0,-500,0,0.001,1,0},
    {0,-500,0,0.0,1,1},
    {0,-INFINITY,0,0.0,1,0}, /* an exactly zero amplitude is still rejected */
    {0,0,0,0.999,1,1},
    {0,-0.5,0,0.1,1,1},
    {0,-0.5,0,0.9,1,0},
    {DBL_MAX,DBL_MAX,-DBL_MAX,0.5,1,1},
    {-DBL_MAX,-DBL_MAX,DBL_MAX,0.5,1,0},
    {NAN,0,0,0.5,0,0},
    {0,NAN,0,0.5,0,0},
    {0,INFINITY,0,0.5,0,0},
    {0,0,-INFINITY,0.5,0,0},
    {0,0,0,NAN,0,0},
    {0,0,0,-0.1,0,0},
    {0,0,0,1.0,0,0},
  };
  for(unsigned i=0;i<sizeof(cases)/sizeof(cases[0]);i++) {
    int accept=-1;
    int valid=BFLogAcceptance(cases[i].p,cases[i].n,cases[i].o,cases[i].u,&accept);
    if(valid!=cases[i].valid || accept!=cases[i].accept) {
      fprintf(stderr,"log acceptance case %u: valid=%d accept=%d\n",i,valid,accept);
      return 1;
    }
  }
  puts("BackFlow log acceptance boundaries passed");
  return 0;
}
