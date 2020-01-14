/*
mVMC - A numerical solver package for a wide range of quantum lattice models based on many-variable Variational Monte Carlo method
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is developed based on the mVMC-mini program
(https://github.com/fiber-miniapp/mVMC-mini)
which follows "The BSD 3-Clause License".

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details. 

You should have received a copy of the GNU General Public License 
along with this program. If not, see http://www.gnu.org/licenses/. 
*/
/*-------------------------------------------------------------
 *[ver.2013.01.21]
 * Unrestircted Hartree-Fock (UHF) method for multi-orbital and
 * multi-interaction Hubbard model
 * Acknowledgements: I thank Daisuke Tahara for his excellent VMC codes.
 * A part of this program is based on his VMC codes.
 * main program
 *-------------------------------------------------------------
 * Copyright (C) 2009- Takahiro MISAWA. All rights reserved.
 *-------------------------------------------------------------*/

#include "Def.h"
#include "matrixlapack.h"
#include "readdef.h"
#include "initial.h"
#include "makeham.h"
#include "diag.h"
#include "green.h"
#include "cal_energy.h"
#include "output.h"
#include "SFMT.h"
#include "xsetmem_def.c"

/*global variables---------------------------------------------*/
struct EDMainCalStruct X;
/*-------------------------------------------------------------*/

int main(int argc, char* argv[]){
    
    /*variable declaration-------------------------------------*/
    //time_t start,mid1,mid2,end;
    char sdt[256];
    FILE *fp;
    int i;
    
    double tmp_eps;

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
    
    if(ReadDefFileNInt(argv[1], &(X.Bind.Def))!=0){
      exit(1);
    };
	
    /*ALLOCATE-------------------------------------------*/
//#include "xsetmem_def.c"
    setmem(&X);
    /*-----------------------------------------------------*/
    
    if(ReadDefFileIdxPara(argv[1], &(X.Bind.Def))!=0){
      exit(1);
    };
    
    /*LARGE VECTORS ARE ALLOCATED*/
#include "xsetmem_large.c"
    /*---------------------------*/
    //Make eps
    tmp_eps=1;
    for(i=0;i<X.Bind.Def.eps_int;i++){
      tmp_eps=tmp_eps*0.1;
    }
    printf("########Input parameters ###########\n");
    printf("tmp_eps=%lf \n",tmp_eps);
    printf("eps_int=%d \n",X.Bind.Def.eps_int);
    printf("mix=%lf \n",X.Bind.Def.mix);
    printf("print=%d \n",X.Bind.Def.print);
    printf("#################################### \n");

    X.Bind.Def.eps=tmp_eps;

	/* initialize Mersenne Twister */
	init_gen_rand(X.Bind.Def.RndSeed);
  //printf("MDEBUG: XXX \n");
	initial(&(X.Bind));
    sprintf(sdt,"%s_check.dat",X.Bind.Def.CDataFileHead);
    fp=fopen(sdt,"w");
    fprintf(fp,"#step,residue,       energy,           # of electrons\n");
    fclose(fp);

    printf("\n########Start: Hartree-Fock calculation ###########\n");
    printf("stp, residue, energy\n");
    for(i=0;i<X.Bind.Def.IterationMax;i++){
      X.Bind.Def.step=i;
     makeham(&(X.Bind));
     diag(&(X.Bind));
     green(&(X.Bind));
     cal_energy(&(X.Bind));
     printf(" %d  %.12lf %.12lf %lf\n",i,X.Bind.Phys.rest,X.Bind.Phys.energy,X.Bind.Phys.num);
     sprintf(sdt,"%s_check.dat",X.Bind.Def.CDataFileHead);
     fp=fopen(sdt,"a");
     fprintf(fp," %d  %.12lf %.12lf %lf\n",i,X.Bind.Phys.rest,X.Bind.Phys.energy,X.Bind.Phys.num);
     fclose(fp);

     if(X.Bind.Phys.rest < X.Bind.Def.eps){
       break;
     } 
    } 

    if(i<X.Bind.Def.IterationMax){
      printf("\nHartree-Fock calculation is finished at %d step. \n\n",i);
    }else{
      printf("\n!! Hartree-Fock calculation is not finished at %d  step!! \n",i);
      printf("\n!! Green functions at the last step are output (these are not converged )!!! \n");
      output(&(X.Bind));
      return -1;
    }
    printf("########Finish: Hartree-Fock calculation ###########\n");

    printf("\n########Start: Calculation of Physical Quantities ###########\n");
    output(&(X.Bind));
    printf("########Finish: Calculation of Physical Quantities ###########\n");


    return 0;
}
