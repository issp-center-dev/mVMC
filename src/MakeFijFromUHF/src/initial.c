#include "initial.h"
void initial(struct BindStruct *X){

    FILE *fp;
    char buf[256],sdt[256];
    int int_i,int_j,int_spin;
    int int_x,int_y;
    int site_i,orb_i;
    int i_max,icnt;
    double sgn,mg,dam,tmp;

    for(int_spin = 0; int_spin < 2; int_spin++){
      for(int_i=0; int_i < X->Def.Nsite; int_i++){
        for(int_j=0; int_j < X->Def.Nsite; int_j++){
          X->Large.G[int_spin][int_i][int_j]   = 0.0;
        }
      }
    }

    for(int_i=0; int_i < X->Def.NInitial; int_i++){
      int_x    = X->Def.Initial[int_i][0]; 
      int_y    = X->Def.Initial[int_i][1]; 
      int_spin = X->Def.Initial[int_i][2]; 
      tmp      = X->Def.ParaInitial[int_i];
      //printf(" Initial: %4d %4d %4d %lf \n",int_x,int_y,int_spin,tmp);
      X->Large.G[int_spin][int_x][int_y]  = tmp ;
    }
}
