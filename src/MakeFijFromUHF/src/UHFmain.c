/*-------------------------------------------------------------
 *[ver.2009.07.08-]
 * Unrestircted Hartree-Fock (UHF) method for multi-orbital and
 * multi-interaction Hubbard model.
 * Acknowledgements: The author thanks Daisuke Tahara for providing him VMC codes.
 * A part of this program is based on Tahara's VMC codes.
 * main program
 *-------------------------------------------------------------
 * Copyright (C) 2009- Takahiro MISAWA. All rights reserved.
 *-------------------------------------------------------------*/

#include "Def.h"
#include "mfmemory.c"
#include "matrixlapack.h"
#include "readdef.h"
#include "check.h"
#include "initial.h"
#include "makeham.h"
#include "diag.h"
#include "green.h"
#include "cal_energy.h"
#include "output.h"
#include "cal_cisajs.h"

double gettimeofday_sec(){
  struct timeval tv;  
  gettimeofday(&tv, NULL);    
  return tv.tv_sec + (double)tv.tv_usec*1e-6;
} 

#define D_FileNameMax 200
/*global variables---------------------------------------------*/
struct EDMainCalStruct X;
/*-------------------------------------------------------------*/

int main(int argc, char* argv[]){
    
    /*variable declaration-------------------------------------*/
    //time_t start,mid1,mid2,end;
    char sdt[256];
	  FILE *fp, *fplist;
	  char configfile[D_FileNameMax];
	  char optfile[D_FileNameMax];
	  char defname[D_FileNameMax];
    
    double t_0,t_1,t_2,t_3,t_4;

	  int mfint[7];/*for malloc*/
	  int i,j,xi;

	  int Nsmp,Nloop;
   	double **vardataBOX, **varADdataBOX;
	  double xtmp,xmax,tmp_eps;
	  int itrstep,jobcont;
    
    //time start
    X.Bind.Time.start=time(NULL);
    
 	  if(argc==1 || argc>3){
		  //ERROR
		  printf("ED Error: *.out(*.exe) NameListFile [OptParaFile]\n");
		  exit(1);
	  }

	  X.Bind.Def.CDataFileHead = (char*)malloc(D_FileNameMax*sizeof(char));
	  X.Bind.Def.CParaFileHead = (char*)malloc(D_FileNameMax*sizeof(char));
	  X.Bind.Def.CPathQtyExe   = (char*)malloc(D_FileNameMax*sizeof(char));
	  X.Bind.Def.CPathAveDev   = (char*)malloc(D_FileNameMax*sizeof(char));
    X.Bind.Def.k_exct = 1;
    X.Bind.Def.nvec   = 1;

	  if(ReadDefFileNInt(argv[1], &(X.Bind.Def))==0){
		  printf("UHF Error: Definition files(*.def) are incomplete.\n");
		  exit(1);
	  };
	  
	  /*ALLOCATE-------------------------------------------*/
	  #include "xsetmem_def.c"
	  /*-----------------------------------------------------*/

	  if(ReadDefFileIdxPara(argv[1], &(X.Bind.Def))==0){
			printf("UHF Error: Indices and Parameters of Definition files(*.def) are incomplete.\n");
		  exit(1);
	  };
    
    check(&(X.Bind)); 
	  /*LARGE VECTORS ARE ALLOCATED*/
	  #include "xsetmem_large.c"
	  /*---------------------------*/
    //Make eps
    tmp_eps=1;
    for(i=0;i<X.Bind.Def.eps_int;i++){
      tmp_eps=tmp_eps*0.1;
    }
    printf("tmp_eps=%lf \n",tmp_eps);
    printf("eps_int=%d \n",X.Bind.Def.eps_int);
    printf("mix=%lf \n",X.Bind.Def.mix);
    printf("print=%d \n",X.Bind.Def.print);
    printf("print=%d \n",X.Bind.Def.print);
    X.Bind.Def.eps=tmp_eps;
    initial(&(X.Bind));
    sprintf(sdt,"%s_check.dat",X.Bind.Def.CDataFileHead);
    fp=fopen(sdt,"w");
    fprintf(fp,"# check start \n");
    fclose(fp);
    for(i=0;i<X.Bind.Def.IterationMax;i++){

     t_0 = gettimeofday_sec();
     makeham(&(X.Bind));
     t_1 = gettimeofday_sec();
     diag(&(X.Bind));
     t_2 = gettimeofday_sec();
     green(&(X.Bind));
     t_3 = gettimeofday_sec();
     calenergy(&(X.Bind));
     t_4 = gettimeofday_sec();
     sprintf(sdt,"%s_check.dat",X.Bind.Def.CDataFileHead);
     fp=fopen(sdt,"a");
     fprintf(fp," %d  %.12lf %.12lf %lf\n",i,X.Bind.Phys.rest,X.Bind.Phys.energy,X.Bind.Phys.num);
     fclose(fp);
     if(X.Bind.Def.print==1){
       printf(" %d  %.12lf %.12lf %lf\n",i,X.Bind.Phys.rest,X.Bind.Phys.energy,X.Bind.Phys.num);
       printf("all: %lf sec\n",t_4-t_0);
       printf("makeham: %lf sec \n",t_1-t_0);
       printf("diag: %lf sec\n",t_2-t_1);
       printf("green: %lf sec\n",t_3-t_2);
       printf("cal: %lf sec\n",t_4-t_3);
     }
      if(X.Bind.Phys.rest < X.Bind.Def.eps){
        break;
      } 
    } 
    output(&(X.Bind));
    cal_cisajs(&(X.Bind));
    if(i<X.Bind.Def.IterationMax){
      printf(" %d finish \n ",i); 
    }else{
      printf(" !! %d Not finish !! \n ",i); 
    } 
    return 0;
}

