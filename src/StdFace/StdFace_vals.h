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
#include <complex.h>

struct StdIntList {
  /*
  Parameters for LATTICE
  */
  double a; /**< The lattice constant */
  double a0;
  double a1;
  int L;
  int W;
  int Lsub;
  int Wsub;
  double Lx;
  double Ly;
  double Wx;
  double Wy;
  int a0L;
  int a0W;
  int a1L;
  int a1W;
  int a0Lsub;
  int a0Wsub;
  int a1Lsub;
  int a1Wsub;
  int S2;
  /*
  Parameters for MODEL
  */
  double mu;
  double complex t;
  double complex tp;
  double complex t0;
  double complex t1;
  double complex t1p;
  double complex t2;
  double complex t2p;
  double U;
  double V;
  double Vp;
  double V0;
  double V1;
  double V1p;
  double V2;
  double V2p;
  /**/
  double JAll;
  double JpAll;
  double J0All;
  double J1All;
  double J1pAll;
  double J2All;
  double J2pAll;
  double J[3][3];
  double Jp[3][3];
  double J0[3][3];
  double J1[3][3];
  double J1p[3][3];
  double J2[3][3];
  double J2p[3][3];
  double D[3][3];
  double h;
  double Gamma;
  double K;
  /*
   Parameters for wavefunction
  */
  int JastrowCut;
  int JastrowAniso;
  /*
   Calculation conditions
  */
  int nelec;
  int NVMCCalMode;
  int NLanczosMode;
  int NDataIdxStart;
  int NDataQtySmp;
  int NSPGaussLeg;
  int NSPStot;
  int NMPTrans;
  int NSROptItrStep;
  int NSROptItrSmp;
  int NSROptFixSmp;
  double DSROptRedCut;
  double DSROptStaDel;
  double DSROptStepDt;
  int NVMCWarmUp;
  int NVMCInterval;
  int NVMCSample;
  int NExUpdatePath;
  int RndSeed;
  int ioutputmode;
  int NSplitSize;
  int NStore;
  int ComplexType;
  /*
   Input strings
  */
  char model[256];
  char lattice[256];
  char outputmode[256];
  char CDataFileHead[256];
  char CParaFileHead[256];
  /*
   Parameter for lattice
  */
  double bW0;
  double bW1;
  double bL0;
  double bL1;
  double bW0sub;
  double bW1sub;
  double bL0sub;
  double bL1sub;
  int NCell;
  int **Cell;
  int NsiteUC;
  double **tau;
  double pi180;
  double phase0;
  double phase1;
  double complex ExpPhase0;
  double complex ExpPhase1;
  int AntiPeriod0;
  int AntiPeriod1;
  /*
   Transfer, Interaction, Locspin
  */
  int nsite;
  int *locspinflag;
  int ntrans;
  int **transindx;
  double complex *trans;
  int nintr;
  int **intrindx;
  double complex *intr;

  int lGC;
  int **Orb;
  int **AntiOrb;
  int NOrb;
  int NSym;
  /*
   Interactions
  */
  int NCintra;
  int **CintraIndx;
  double *Cintra;
  int NCinter;
  int **CinterIndx;
  double *Cinter;
  int NHund;
  int **HundIndx;
  double *Hund;
  int NEx;
  int **ExIndx;
  double *Ex;
};
