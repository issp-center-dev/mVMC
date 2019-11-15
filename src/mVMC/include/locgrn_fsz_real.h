#ifndef _LOCGRN_FSZ_REAL
#define _LOCGRN_FSZ_REAL

double GreenFunc1_fsz_real(const int ri, const int rj, const int s, const double ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,int *eleSpn,
                  int *projCntNew, double *buffer);

double GreenFunc1_fsz2_real(const int ri, const int rj, const int s,const int t, const double ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,int *eleSpn,
                  int *projCntNew, double *buffer);

double GreenFunc2_fsz_real(const int ri, const int rj, const int rk, const int rl,
                  const int s, const int t, const double ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,int *eleSpn,
                  int *projCntNew, double *buffer);


#endif
