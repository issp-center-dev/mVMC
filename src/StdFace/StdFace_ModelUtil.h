/*
HPhi-mVMC-StdFace - Common input generator
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
#include <complex.h>
#include <stdio.h>

void StdFace_exit(int errorcode);

void StdFace_intr(struct StdIntList *StdI, double complex intr0,
  int site1, int spin1, int site2, int spin2,
  int site3, int spin3, int site4, int spin4);

void StdFace_Hopping(struct StdIntList *StdI, double complex trans0, int isite, int jsite, double *dR);
void StdFace_trans(struct StdIntList *StdI,double complex trans0,int isite,int ispin,int jsite,int jspin);
void StdFace_HubbardLocal(struct StdIntList *StdI, double mu0, double h0,
  double Gamma0, double U0, int isite);
void StdFace_MagField(struct StdIntList *StdI, int S2, double h, double Gamma, int isite);

void StdFace_Coulomb(struct StdIntList *StdI, double V, int isite, int jsite);
void StdFace_GeneralJ(struct StdIntList *StdI, double J[3][3],
  int Si2, int Sj2, int isite, int jsite);

void StdFace_PrintVal_d(char* valname, double *val, double val0);
void StdFace_PrintVal_dd(char* valname, double *val, double val0, double val1);
void StdFace_PrintVal_c(char* valname, double complex *val, double complex val0);
void StdFace_PrintVal_i(char* valname, int *val, int val0);

void StdFace_NotUsed_d(char* valname, double val);
void StdFace_NotUsed_i(char* valname, int val);
void StdFace_NotUsed_c(char* valname, double complex val);
void StdFace_NotUsed_J(char* valname, double JAll, double J[3][3]);

void StdFace_RequiredVal_i(char* valname, int val);
void StdFace_InputSpinNN(double J[3][3], double JAll, double J0[3][3], double J0All, char *J0name);
void StdFace_InputSpin(double Jp[3][3], double JpAll, char *Jpname);
void StdFace_InputCoulombV(double V, double *V0, char *V0name);
void StdFace_InputHopp(double complex t, double complex *t0, char *t0name);

void StdFace_InitSite(struct StdIntList *StdI, FILE *fp, int dim);
void StdFace_SetLabel(struct StdIntList *StdI, FILE *fp,
  int iW, int iL, int diW, int diL, int isiteUC, int jsiteUC,
  int *isite, int *jsite, int connect, double complex *Cphase, double *dR);
void StdFace_PrintGeometry(struct StdIntList *StdI);
void StdFace_MallocInteractions(struct StdIntList *StdI, int ntransMax, int nintrMax);
void StdFace_FindSite(struct StdIntList *StdI,
  int iW, int iL, int iH, int diW, int diL, int diH,
  int isiteUC, int jsiteUC,
  int *isite, int *jsite, double complex *Cphase, double *dR);
void StdFace_PrintXSF(struct StdIntList *StdI);

void StdFace_Tetragonal(struct StdIntList *StdI);
void StdFace_Chain(struct StdIntList *StdI);
void StdFace_Ladder(struct StdIntList *StdI);
void StdFace_Triangular(struct StdIntList *StdI);
void StdFace_Honeycomb(struct StdIntList *StdI);
void StdFace_Kagome(struct StdIntList *StdI);
void StdFace_Orthorhombic(struct StdIntList *StdI);
void StdFace_FCOrtho(struct StdIntList *StdI);
void StdFace_Pyrochlore(struct StdIntList *StdI);
void StdFace_Wannier90(struct StdIntList *StdI);

#if defined(_HPhi)
void StdFace_Chain_Boost(struct StdIntList *StdI);
void StdFace_Ladder_Boost(struct StdIntList *StdI);
void StdFace_Honeycomb_Boost(struct StdIntList *StdI);
void StdFace_Kagome_Boost(struct StdIntList *StdI);
#elif defined(_mVMC)
void StdFace_generate_orb(struct StdIntList *StdI);
void StdFace_Proj(struct StdIntList *StdI);
void PrintJastrow(struct StdIntList *StdI);
#endif
