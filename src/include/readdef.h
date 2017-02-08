/*
mVMC - A numerical solver package for a wide range of quantum lattice models based on many-variable Variational Monte Carlo method
Copyright (C) 2016 Takahiro Misawa, Satoshi Morita, Takahiro Ohgoe, Kota Ido, Mitsuaki Kawamura, Takeo Kato, Masatoshi Imada.

his program is developed based on the mVMC-mini program
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
#pragma once
#define D_FileNameMaxReadDef 256 /*!<  Max length of words for file name*/
#define D_CharTmpReadDef     200 /*!<  Max length of reading words from input files*/
#define D_CharKWDMAX     200 /*!<  Max length of words for keyword*/

/*!< Number of ignore lines in def files */
#define IgnoreLinesInDef 5

/**
 * Keyword List in NameListFile.
 **/
static char cKWListOfFileNameList[][D_CharTmpReadDef]={
  "ModPara", "LocSpin",
  "Trans", "CoulombIntra", "CoulombInter",
  "Hund", "PairHop", "Exchange",
  "Gutzwiller", "Jastrow",
  "DH2", "DH4", "Orbital",
  "TransSym", "InGutzwiller", "InJastrow",
  "InDH2", "InDH4", "InOrbital",
  "OneBodyG", "TwoBodyG", "TwoBodyGEx",
  "InterAll", "OptTrans", "InOptTrans"
};

/**
 * Number of Keyword List in NameListFile for this program.  
 **/
enum KWIdxInt{
  KWModPara, KWLocSpin,
  KWTrans, KWCoulombIntra,KWCoulombInter,
  KWHund, KWPairHop, KWExchange,
  KWGutzwiller, KWJastrow,
  KWDH2, KWDH4, KWOrbital,
  KWTransSym, KWInGutzwiller, KWInJastrow,
  KWInDH2, KWInDH4, KWInOrbital,
  KWOneBodyG, KWTwoBodyG, KWTwoBodyGEx,
  KWInterAll, KWOptTrans, KWInOptTrans,
  KWIdxInt_end
};

/**
 * File Name List in NameListFile.
 **/
static char (*cFileNameListFile)[D_CharTmpReadDef];


enum ParamIdxInt{
  IdxVMCCalcMode, IdxLanczosMode, IdxDataIdxStart, 
  IdxDataQtySmp, IdxNsite, IdxNe, 
  IdxSPGaussLeg, IdxSPStot, IdxMPTrans,
  IdxSROptItrStep, IdxSROptItrSmp, IdxSROptFixSmp,
  IdxVMCWarmUp, IdxVMCInterval, IdxVMCSample,
  IdxExUpdatePath, IdxRndSeed, IdxSplitSize,
  IdxNLocSpin,IdxNTrans,IdxNCoulombIntra,
  IdxNCoulombInter, IdxNHund, IdxNPairHop, 
  IdxNExchange, IdxNGutz, IdxNJast,
  IdxNDH2, IdxNDH4, IdxNOrbit,
  IdxNQPTrans, IdxNOneBodyG, IdxNTwoBodyG,
  IdxNTwoBodyGEx, IdxNInterAll, IdxNQPOptTrans,
  IdxSROptCGMaxIter,
  ParamIdxInt_End
};

enum ParamIdxDouble{
  IdxSROptRedCut, IdxSROptStaDel, IdxSROptStepDt,
  IdxSROptCGTol,
  ParamIdxDouble_End
};


int CheckWords( const char* ctmp, const char* cKeyWord);
int CheckKW(const char* cKW, char  cKWList[][D_CharTmpReadDef], int iSizeOfKWidx, int* iKWidx);
int GetKWWithIdx(char *ctmpLine, char *ctmp, int *itmp);
int ValidateValue(const int icheckValue, const int ilowestValue, const int iHighestValue);
int GetFileName(const char* cFileListNameFile, char cFileNameList[][D_CharTmpReadDef]);

void SetDefultValuesModPara(int *buf, double* bufDouble);

int ReadInputParameters(char *xNameListFile, MPI_Comm comm);

// Flag for giving type of optimized parameters
// 0: real, 1: complex
int iComplexFlgGutzwiller=0;
int iComplexFlgJastrow=0;
int iComplexFlgDH2=0;
int iComplexFlgDH4=0;
int iComplexFlgOrbital=0;
