/*-------------------------------------------------------------
 *[ver.2009.3.31]
 * Exact Diagonalization with Lanczos Method
 *-------------------------------------------------------------
 * Copyright (C) 2006-2009 Takahiro MISAWA. All rights reserved.
 *-------------------------------------------------------------*/

/*=================================================================================================*/
struct DefineList{
	char *CDataFileHead, *CParaFileHead, *CPathQtyExe, *CPathAveDev;
  int nvec,k_exct;
	int NDataIdxStart, NDataQtySmp;

	int NVMCCalMode, NLanczosMode;

	int Nsite, Ne, Nsize, IterationMax;
  int print;
  int eps_int;
  double eps,mix;
  int St;
	
	int *LocSpn, NLocSpn, **All_pair;

	int *OptFlag, fidx;
  unsigned int *Tpow;

  int *EDChemi, EDNChemi; double *EDParaChemi;
	int **Transfer, **EDTransfer, NTransfer,EDNTransfer;					double *ParaTransfer,*EDParaTransfer;
	int **CoulombIntra, NCoulombIntra;			double *ParaCoulombIntra;
	int **CoulombInter, NCoulombInter;			double *ParaCoulombInter;
	int **HundCoupling, NHundCoupling;			double *ParaHundCoupling;
	int **PairHopping, NPairHopping;			double *ParaPairHopping;
	int **ExchangeCoupling, NExchangeCoupling;	double *ParaExchangeCoupling;
  int **Initial,NInitial;  
  double *ParaInitial;
  

  int **CisAjs, NCisAjs;
  int **CisAjsCktAltDC, NCisAjsCktAltDC;

	
};

struct CheckList{
  int idim_max,sdim; 
  double max_mem;
};

struct LargeList{
    double ***Ham,***G,***G_old; 
    double **R_SLT_0,**L_SLT_0;
    double **R_SLT_1,**L_SLT_1;
    double **EigenValues;
    double *tmp;
};

struct PhysList{
	double energy,doublon;
  double rest,num;
	double *spin_real_cor;
  double *charge_real_cor;
  double *loc_spin_z;
  double Target_energy;
};	
struct TimeList{
  time_t start,mid1,mid2,end;
};
/*=================================================================================================*/
struct BindStruct{
	struct DefineList  Def;
	struct CheckList   Check;
	struct LargeList   Large;
	struct PhysList    Phys;
	struct TimeList    Time;
};
/*=================================================================================================*/
struct EDMainCalStruct{
	struct BindStruct Bind;
};


