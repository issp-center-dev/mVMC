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
 *[ver.2009.3.31]
 * Exact Diagonalization with Lanczos Method
 *-------------------------------------------------------------
 * Copyright (C) 2006-2009 Takahiro MISAWA. All rights reserved.
 *-------------------------------------------------------------*/

/*=================================================================================================*/
struct DefineList {
    int step;
    char *CDataFileHead, *CParaFileHead, *CPathQtyExe, *CPathAveDev;
    int nvec, k_exct;
    int NDataIdxStart, NDataQtySmp;

    int NVMCCalMode, NLanczosMode;

    int Nsite, Ne, Nsize, IterationMax;
    int TwoSz, Ncond;
    int print;
    int eps_int;
    double eps, mix;
    int St;
    int RndSeed;
    int NMPTrans;
    int eps_int_slater;
    int APFlag;
    int *LocSpn, NLocSpn, **All_pair;

    int *OptFlag, fidx;
    unsigned int *Tpow;

    int *EDChemi, EDNChemi;
    double *EDParaChemi;
    int **Transfer, **EDTransfer, NTransfer, EDNTransfer;
    complex double *ParaTransfer;
    int **CoulombIntra, NCoulombIntra;
    double *ParaCoulombIntra;
    int **CoulombInter, NCoulombInter;
    double *ParaCoulombInter;
    int **HundCoupling, NHundCoupling;
    double *ParaHundCoupling;
    int **PairHopping, NPairHopping;
    double *ParaPairHopping;
    int **ExchangeCoupling, NExchangeCoupling;
    double *ParaExchangeCoupling;
    int NOrbitalIdx, **OrbitalIdx; /* [Nsite][Nsite] */
    int **OrbitalSgn; /* OrbitalSgn[2*Nsite][2*Nsite] = +1 or -1 */
    int **Initial, NInitial;
    complex double *ParaInitial;
    //double *ParaInitial;
    //double *ParaInitial_theta;
    int NOrbitalAP, NOrbitalP;


    int **CisAjs, NCisAjs;
    int **CisAjsCktAltDC, NCisAjsCktAltDC;

    int iFlgOrbitalGeneral;
    int OrbitalOutputMode;


};

struct CheckList {
    int idim_max, sdim;
    double max_mem;
};

struct LargeList {
    double complex **Ham, **G, **G_old;
    double complex **R_SLT, **L_SLT;
    double *EigenValues;
    double *tmp;
};

struct PhysList {
    double energy, doublon;
    double rest, num;
    double *spin_real_cor;
    double *charge_real_cor;
    double *loc_spin_z;
    double Target_energy;
};

struct TimeList {
    time_t start, mid1, mid2, end;
};

/*=================================================================================================*/
struct BindStruct {
    struct DefineList Def;
    struct CheckList Check;
    struct LargeList Large;
    struct PhysList Phys;
    struct TimeList Time;
};

/*=================================================================================================*/
struct EDMainCalStruct {
    struct BindStruct Bind;
};


