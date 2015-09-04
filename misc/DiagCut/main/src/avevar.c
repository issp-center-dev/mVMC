/*-------------------------------------------------------------
 * Variational Monte Carlo
 * calculate average and variance of variational parameters
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/
void StoreOptData(int sample);
void OutputOptData();

void StoreOptData(int sample){
  const int n = 2+NPara;
  int i;
  double *optData = SROptData + n*sample;

  optData[0] = Etot;
  optData[1] = Etot2;
  for(i=0;i<NPara;i++) optData[i+2] = Para[i];
  
  return;
}

void OutputOptData() {
  const int n = 2+NPara;
  char fileName[D_FileNameMax];
  FILE *fp;
  double ave,var,data;
  int sample,i;

  sprintf(fileName, "%s_opt.dat", CParaFileHead);
  fp = fopen(fileName, "w");

  if(NSROptItrSmp==1) {
    for(i=0;i<n;i++) {
      fprintf(fp,"% .18e % .18e ", SROptData[i], 0.0);
    }
  } else {
    for(i=0;i<n;i++) {
      ave=0.0;
      for(sample=0;sample<NSROptItrSmp;sample++) {
        ave += SROptData[i+n*sample];
      }
      ave /= (double)(NSROptItrSmp);

      var = 0.0;
      for(sample=0;sample<NSROptItrSmp;sample++) {
        data = SROptData[i+n*sample] - ave;
        var += data*data;
      }
      var = sqrt( var/((double)(NSROptItrSmp)-1.0) );

      fprintf(fp,"% .18e % .18e ", ave, var);
    }
  }
  fprintf(fp, "\n");
  fclose(fp);

  return;
}
