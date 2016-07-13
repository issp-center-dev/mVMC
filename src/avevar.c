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

void OutputGutzwiller(FILE *fp_all, int count_i){
  FILE *fp_gutz;
  fp_gutz = fopen("gutzwiller_opt.out", "w");
  WriteHeader("NGutzwillerIdx", NGutzwillerIdx,  fp_gutz);

  for(i=0; i<NGutzwillerIdx; i++){
    CalcAveVar(count_i+i, n, &ave, &var);  
    fprintf(fp_gutz, "%d % .18e % .18e \n",
            GutzwillerIdx[i], creal(ave), cimag(ave));
    fprintf(fp_all,"% .18e % .18e % .18e ",
            creal(ave), cimag(ave), var);
  }
  fp_gutz.close();
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
    OutputGutzwiller(fp, count_i);

    count_i += NGutzwillerIdx;
    for(i=0; i<NJastrowIdx; i++){
      CalcAveVar(count_i+i, n, &ave, &var);
      fprintf(fp,"% .18e % .18e % .18e ",
              creal(ave), cimag(ave), var);
    }

    count_i += NJastrowIdx;
    for(i=0; i<2*3*NDoublonHolon2siteIdx; i++){
      CalcAveVar(count_i+i, n, &ave, &var);
      fprintf(fp,"% .18e % .18e % .18e ",
              creal(ave), cimag(ave), var);
    }

    count_i +=2*3*NDoublonHolon2siteIdx;
    for(i=0; i<2*5*NDoublonHolon4siteIdx; i++){
      CalcAveVar(count_i+i, n, &ave, &var);
      fprintf(fp,"% .18e % .18e % .18e ",
              creal(ave), cimag(ave), var);
    }

    count_i +=2*5*NDoublonHolon4siteIdx;
    for(i=0; i<NSlater; i++){
      CalcAveVar(count_i+i, n, &ave, &var);
      fprintf(fp,"% .18e % .18e % .18e ",
              creal(ave), cimag(ave), var);
    }

    count_i +=NSlater;
    for(i=0; i<NOptTrans; i++){
      CalcAveVar(count_i+i, n, &ave, &var);
      fprintf(fp,"% .18e % .18e % .18e ",
              creal(ave), cimag(ave), var);
    }

    count_i += NOptTrans;
    if(count_i != n){
      printf("Debug: NProj=%d, NPara=%d\n", NProj, NPara);
      printf("Debug: Error %d, %d\n", count_i, n);
    }    
  }
  fprintf(fp, "\n");
  fclose(fp);

  return;
}

void OutputJastrow(){
  
}

void OutputDoublonHolon2site(){
  
}

void OutputDoublonHolon4site(){
  
}

void OutputSlater(){
}

void OutputOptTrans(){
}
