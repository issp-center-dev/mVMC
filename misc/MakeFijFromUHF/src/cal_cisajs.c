void cal_cisajs(struct BindStruct *X){

    time_t start,end;
    FILE *fp;
    char sdt[256];
    int int_i,int_spin,site_1,site_2;
    double tmp;
 
    sprintf(sdt, "%s_UHF_cisajs.dat", X->Def.CDataFileHead);
    fp = fopen(sdt,"w");

    for(int_i=0; int_i< X->Def.NCisAjs; int_i++){
            
      site_1    = X->Def.CisAjs[int_i][0];
      site_2    = X->Def.CisAjs[int_i][1];
      int_spin  = X->Def.CisAjs[int_i][2];
      tmp       = X->Large.G[int_spin][site_1][site_2];

      fprintf(fp," %4d %4d %4d %.10lf \n",site_1,site_2,int_spin,tmp);
    }
    fclose(fp);
}
