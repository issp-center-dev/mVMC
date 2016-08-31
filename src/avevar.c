/*-------------------------------------------------------------
 * Variational Monte Carlo
 * calculate average and variance of variational parameters
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/
#include "./include/readdef.h"

void StoreOptData(int sample);
void OutputOptData();

void CalcAveVar(int i, int n, double complex* _ave, double* _var){
  int sample;
  double complex data;
  double complex ave=0;
  double var=0;
  
  for(sample=0;sample<NSROptItrSmp;sample++) {
    ave += SROptData[i+n*sample];
  }
  ave /= (double)(NSROptItrSmp);
  
  var = 0.0;
  for(sample=0;sample<NSROptItrSmp;sample++) {
    data = SROptData[i+n*sample] - ave;
    var += creal(data*conj(data));//TBC
  }
  var = sqrt( var/((double)(NSROptItrSmp)-1.0) );

  *_ave= ave;
  *_var= var;
}

void WriteHeader(char* cNKWidx, int NKWidx, FILE *fp){
  fprintf(fp, "======================\n");
  fprintf(fp, cNKWidx);
  fprintf(fp, "  %d\n", NKWidx);
  fprintf(fp, "======================\n");
  fprintf(fp, "======================\n");
  fprintf(fp, "======================\n");    
}

void Child_OutputOptData(FILE *fp_all, char* cFileName, char* cNKWidx, int Nidx_head, int Nidx, int count_i, int n){
  FILE *fp_out;
  int i;
  double complex ave;
  double var;
  fp_out = fopen(cFileName, "w");
  WriteHeader(cNKWidx, Nidx_head,  fp_out);
  for(i=0; i<Nidx; i++){
    CalcAveVar(count_i+i, n, &ave, &var);  
    fprintf(fp_out, "%d % .18e % .18e \n",
            i, creal(ave), cimag(ave));
    fprintf(fp_all,"% .18e % .18e % .18e ",
            creal(ave), cimag(ave), var);
  }
  fclose(fp_out);
}
  
void StoreOptData(int sample){
  const int n = 2+NPara;
  int i;
  double complex *optData = SROptData + n*sample;

  optData[0] = Etot;
  optData[1] = Etot2;
  for(i=0;i<NPara;i++) optData[i+2] = Para[i];
  
  return;
}

void OutputOptData() {
  const int n = 2+NPara;
  char fileName[D_FileNameMax];
  FILE *fp;
  //double ave,var,data;
  double complex ave,data;
  double var;
  int sample,i;
  int count_i;
  sprintf(fileName, "%s_opt.dat", CParaFileHead);
  fp = fopen(fileName, "w");

  if(NSROptItrSmp==1) {
    for(i=0;i<n;i++) {
      fprintf(fp,"% .18e % .18e ", creal(SROptData[i]), 0.0);//TBC
    }
  } else {    
    //output <H> and <H^2>
    for(i=0;i<2;i++) {
      CalcAveVar(i, n, &ave, &var);
      fprintf(fp,"% .18e % .18e % .18e ",
              creal(ave), cimag(ave), var);
    }

    count_i=2;
    if(NGutzwillerIdx !=0){
      sprintf(fileName, "%s_gutzwiller_opt.dat", CParaFileHead);
      Child_OutputOptData(fp, fileName, "NGutzwillerIdx",
			  NGutzwillerIdx, NGutzwillerIdx, count_i, n);
      count_i += NGutzwillerIdx;
    }

    if(NJastrowIdx !=0){
      sprintf(fileName, "%s_jastrow_opt.dat", CParaFileHead);
      Child_OutputOptData(fp, fileName, "NJastrowIdx",
			  NJastrowIdx, NJastrowIdx, count_i, n);
      count_i += NJastrowIdx;
    }

    if(NDoublonHolon2siteIdx != 0){
      sprintf(fileName, "%s_doublonHolon2site_opt.dat", CParaFileHead);
      Child_OutputOptData(fp, fileName, "NDoublonHolon2siteIdx",
			  NDoublonHolon2siteIdx, NDoublonHolon2siteIdx*2*3, 
			  count_i, n);
      count_i +=2*3*NDoublonHolon2siteIdx;
    }

    if(NDoublonHolon4siteIdx != 0){
      sprintf(fileName, "%s_doublonHolon4site_opt.dat", CParaFileHead);
      Child_OutputOptData(fp, fileName, "NDoublonHolon4siteIdx",
			  NDoublonHolon4siteIdx, NDoublonHolon4siteIdx*2*5, 
			  count_i, n);
      count_i +=2*5*NDoublonHolon4siteIdx;
    }

    if(NSlater != 0){
      sprintf(fileName, "%s_orbital_opt.dat", CParaFileHead);
      Child_OutputOptData(fp, fileName, "NOrbitalIdx",
			  NSlater, NSlater, count_i, n);
      count_i +=NSlater;
    }

    if(NOptTrans !=0){
      sprintf(fileName, "%s_trans_opt.dat", CParaFileHead);
      Child_OutputOptData(fp, fileName, "NQPOptTrans",
			  NOptTrans, NOptTrans, count_i, n);
      count_i += NOptTrans;
    }
  }
  fprintf(fp, "\n");
  fclose(fp);

  return;
}
