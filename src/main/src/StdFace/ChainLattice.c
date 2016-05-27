/*
HPhi  -  Quantum Lattice Model Simulator
Copyright (C) 2015 The University of Tokyo

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "StdFace_vals.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "StdFace_ModelUtil.h"
#include <complex.h>
#include "../include/wrapperMPI.h"
#include <string.h>

/**
 *
 * Setup a Hamiltonian for the Hubbard model on a Chain lattice
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void StdFace_Chain(struct StdIntList *StdI, char *model)
{
  FILE *fp;
  int isite, jsite;
  int iL;
  int ktrans, kintr;
  /**/
  fprintf(stdout, "\n");
  fprintf(stdout, "#######  Parameter Summary  #######\n");
  fprintf(stdout, "\n");
  /**/
  StdI->NsiteUC = 1;
  fprintf(stdout, "  @ Lattice Size & Shape\n\n");
  StdFace_RequiredVal_i("L", StdI->L);
  StdFace_NotUsed_i("W", StdI->W);
  StdI->W = 1;
  StdFace_NotUsed_i("a0W", StdI->a0W);
  StdFace_NotUsed_i("a0L", StdI->a0L);
  StdFace_NotUsed_i("a1W", StdI->a1W);
  StdFace_NotUsed_i("a1L", StdI->a1L);
  /**/
  StdI->a0 = 1.0; StdI->a1 = 1.0;
  fp = fopen("/dev/null", "w");
  StdFace_InitSite2D(StdI, fp, StdI->a0, 0.0, 0.0, StdI->a1);
  fclose(fp);
  StdI->tau[0][0] = 0.0; StdI->tau[0][1] = 0.0;
  /**/
  fprintf(stdout, "\n  @ Hamiltonian \n\n");
  StdFace_NotUsed_J("J1", StdI->J1All, StdI->J1);
  StdFace_NotUsed_J("J2", StdI->J2All, StdI->J2);
  StdFace_NotUsed_J("J1'", StdI->J1pAll, StdI->J1p);
  StdFace_NotUsed_J("J2'", StdI->J2pAll, StdI->J2p);
  StdFace_NotUsed_c("t1", StdI->t1);
  StdFace_NotUsed_c("t2", StdI->t2);
  StdFace_NotUsed_d("t1'", StdI->t1p);
  StdFace_NotUsed_d("t2'", StdI->t2p);
  StdFace_NotUsed_d("V1", StdI->V1);
  StdFace_NotUsed_d("V2", StdI->V2);
  StdFace_NotUsed_d("V1'", StdI->V1p);
  StdFace_NotUsed_d("V2'", StdI->V2p);
  StdFace_NotUsed_d("K", StdI->K);
  /**/
  if (strcmp(StdI->model, "spin") == 0 ) {
    StdFace_PrintVal_i("2S", &StdI->S2, 1);
    StdFace_PrintVal_d("h", &StdI->h, 0.0);
    StdFace_PrintVal_d("Gamma", &StdI->Gamma, 0.0);
    StdFace_PrintVal_d("D", &StdI->D[2][2], 0.0);
    StdFace_InputSpinNN(StdI, StdI->J0, StdI->J0All, "J0");
    StdFace_InputSpin(StdI, StdI->Jp, StdI->JpAll, "J'");
    /**/
    StdFace_NotUsed_d("mu", StdI->mu);
    StdFace_NotUsed_d("U", StdI->U);
    StdFace_NotUsed_c("t", StdI->t);
    StdFace_NotUsed_c("t0", StdI->t0);
    StdFace_NotUsed_c("t'", StdI->tp);
    StdFace_NotUsed_d("V", StdI->V);
    StdFace_NotUsed_d("V0", StdI->V0);
    StdFace_NotUsed_d("V'", StdI->Vp);
  }/*if (strcmp(StdI->model, "spin") == 0 )*/
  else {
    StdFace_PrintVal_d("mu", &StdI->mu, 0.0);
    StdFace_PrintVal_d("U", &StdI->U, 0.0);
    StdFace_InputHopp(StdI, &StdI->t0, "t0");
    StdFace_PrintVal_c("t'", &StdI->tp, 0.0);
    StdFace_InputCoulombV(StdI, &StdI->V0, "V0");
    StdFace_PrintVal_d("V'", &StdI->Vp, 0.0);

    StdFace_NotUsed_J("J0", StdI->J0All, StdI->J0);
    StdFace_NotUsed_J("J'", StdI->JpAll, StdI->Jp);
    StdFace_NotUsed_d("h", StdI->h);
    StdFace_NotUsed_d("Gamma", StdI->Gamma);
    StdFace_NotUsed_d("D", StdI->D[2][2]);

    if (strcmp(StdI->model, "hubbard") == 0 ) {
      StdFace_NotUsed_i("2S", StdI->S2);
      StdFace_NotUsed_J("J", StdI->JAll, StdI->J);
    }
    else if (strcmp(StdI->model, "kondo") == 0 ) {
      StdFace_PrintVal_i("2S", &StdI->S2, 1);
      StdFace_InputSpin(StdI, StdI->J, StdI->JAll, "J");
    }
  }/*if (strcmp(StdI->model, "spin") != 0 )*/
  fprintf(stdout, "\n  @ Numerical conditions\n\n");
  /*
  Local Spin
  */
  StdI->nsite = StdI->L;
  if (strcmp(StdI->model, "kondo") == 0 ) StdI->nsite *= 2;
  StdI->locspinflag = (int *)malloc(sizeof(int) * StdI->nsite);
  /**/
  if (strcmp(StdI->model, "spin") == 0 )
    for (isite = 0; isite < StdI->nsite; isite++)StdI->locspinflag[isite] = StdI->S2;
  else if (strcmp(StdI->model, "hubbard") == 0 ) 
    for (isite = 0; isite < StdI->nsite; isite++)StdI->locspinflag[isite] = 0;
  else if (strcmp(StdI->model, "kondo") == 0 ) 
    for (isite = 0; isite < StdI->nsite / 2; isite++) {
      StdI->locspinflag[isite] = StdI->S2;
      StdI->locspinflag[isite + StdI->nsite / 2] = 0;
    }
  /*
   The number of Transfer & Interaction
  */
  if (strcmp(StdI->model, "spin") == 0 ) {
    StdI->ntrans = StdI->L * (StdI->S2 + 1/*h*/ + 2 * StdI->S2/*Gamma*/);
    StdI->nintr = StdI->L * (StdI->NsiteUC/*D*/ + 1/*J*/ + 1/*J'*/)
      * (3 * StdI->S2 + 1) * (3 * StdI->S2 + 1);
  }
  else {
    StdI->ntrans = StdI->L * 2/*spin*/ * (StdI->NsiteUC/*mu*/ + 2/*t*/ + 2/*t'*/);
    StdI->nintr = StdI->L * (StdI->NsiteUC/*U*/ + 4 * (1/*V*/ + 1/*V'*/));

    if(strcmp(StdI->model, "kondo") == 0 ) 
      StdI->nintr += StdI->nsite / 2 * (3 * 1 + 1) * (3 * StdI->S2 + 1);
  }
  /**/
  StdI->transindx = (int **)malloc(sizeof(int*) * StdI->ntrans);
  StdI->trans = (double complex *)malloc(sizeof(double complex) * StdI->ntrans);
  for (ktrans = 0; ktrans < StdI->ntrans; ktrans++){
    StdI->transindx[ktrans] = (int *)malloc(sizeof(int) * 4);
  }
  /**/
  StdI->intrindx = (int **)malloc(sizeof(int*) * StdI->nintr);
  StdI->intr = (double complex *)malloc(sizeof(double complex) * StdI->nintr);
  for (kintr = 0; kintr < StdI->nintr; kintr++) {
    StdI->intrindx[kintr] = (int *)malloc(sizeof(int) * 8);
  }
  /*
   Set Transfer & Interaction
  */
  StdI->ntrans = 0;
  StdI->nintr = 0;
  for (iL = 0; iL < StdI->L; iL++){

    isite = iL;
    if (strcmp(StdI->model, "kondo") == 0 ) isite += StdI->L;
    /*
     Local term
    */
    if (strcmp(StdI->model, "spin") == 0 ) {
      StdFace_MagField(StdI, StdI->S2, -StdI->h, -StdI->Gamma, isite);
      StdFace_GeneralJ(StdI, StdI->D, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, StdI->mu, isite, isite);
      StdFace_intr(StdI, StdI->U, isite, 0, isite, 0, isite, 1, isite, 1);
      /**/
      if (strcmp(StdI->model, "kondo") == 0 ) {
        jsite = iL;
        StdFace_GeneralJ(StdI, StdI->J, 1, StdI->S2, isite, jsite);
      }/*if (strcmp(StdI->model, "kondo") == 0 )*/
    }/*if (model != "spin")*/
    /*
    Nearest neighbor
   */
    jsite = (iL + 1) % StdI->L;
    if (strcmp(StdI->model, "kondo") == 0 ) jsite += StdI->L;
    /**/
    if (strcmp(StdI->model, "spin") == 0 ) {
      StdFace_GeneralJ(StdI, StdI->J0, StdI->S2, StdI->S2, isite, jsite);
    }
    else {
      StdFace_Hopping(StdI, StdI->t0, isite, jsite);
      StdFace_Coulomb(StdI, StdI->V0, isite, jsite);
    }
    /*
    Second nearest neighbor
    */
    jsite = (iL + 2) % StdI->L;
    if (strcmp(StdI->model, "kondo") == 0 ) jsite += StdI->L;
    /**/
    if (strcmp(StdI->model, "spin") == 0 ) {
      StdFace_GeneralJ(StdI, StdI->Jp, StdI->S2, StdI->S2, isite, jsite);
    }
    else {
      StdFace_Hopping(StdI, StdI->tp, isite, jsite);
      StdFace_Coulomb(StdI, StdI->Vp, isite, jsite);
    }
  }/*for (iL = 0; iL < StdI->L; iL++)*/
   /*
   Set Orbital index
   */
  generate_orb(StdI);
}
