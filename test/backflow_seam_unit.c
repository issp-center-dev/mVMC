/* Strict BFRange parsing, packed allocation and active transformation contract. */
#include <sys/stat.h>
#include <unistd.h>
#include "../src/mVMC/backflow.c"

#define CHECK(x) do { if (!(x)) { fprintf(stderr, "FAIL line %d: %s\n", __LINE__, #x); exit(1); } } while (0)

static FILE *input(const char *body, size_t size) {
  FILE *fp = fopen("tmp/backflow_seam_unit.input", "w+b");
  CHECK(fp != NULL);
  CHECK(fwrite(body, 1, size, fp) == size);
  rewind(fp);
  return fp;
}

static int read_range(int columns, int ap, int badrow, const char *replacement,
                      const char *trailer) {
  char body[4096];
  size_t used = 0;
  int i, r, row = 0, result;
  FILE *fp;
  APFlag = ap;
  for (i = 0; i < 5; i++) used += snprintf(body + used, sizeof(body)-used, "header\n");
  for (i = 0; i < 4; i++) {
    for (r = 0; r < 3; r++, row++) {
      const int k = (i + (r == 2 ? 3 : r)) % 4;
      if (row == badrow) {
        used += snprintf(body + used, sizeof(body)-used, "%s", replacement);
      } else {
        used += snprintf(body + used, sizeof(body)-used, "%d %d %d", i, k, r != 0);
        if (columns == 4)
          used += snprintf(body + used, sizeof(body)-used, " %d", ap && abs(i-k)==3 ? -1 : 1);
        used += snprintf(body + used, sizeof(body)-used, "\n");
      }
    }
  }
  if (trailer) used += snprintf(body+used, sizeof(body)-used, "%s", trailer);
  CHECK(used < sizeof(body));
  fp = input(body, used);
  result = BFReadRange(fp, "unit");
  fclose(fp);
  return result;
}

int main(void) {
  int storage[46], expected[44], *cursor = storage + 1;
  int t0[4]={0,1,2,3}, t1[4]={1,2,3,0};
  int s0[4]={1,1,1,1}, s1[4]={1,1,1,-1};
  int *trans[2]={t0,t1}, *signs[2]={s0,s1};
  int *opt[1]={t1}, *optsigns[1]={s1};
  const char *bad[] = {"0 1\n", "0 1 1 1 1\n", "0 1 1\n", "0 1 1 0\n",
    "0 1 1 2\n", "0 1 1 1.0\n", "0 1 1 1junk\n", "0 1 1 2147483648\n",
    "0 1 1 -2147483649\n", "\n", "# comment\n", "0 1 1 1 #comment\n",
    "-1 1 1 1\n", "4 1 1 1\n", "0 1 -1 1\n", "0 1 2 1\n"};
  const int created = mkdir("tmp", 0700) == 0;
  size_t i;
  Nsite=4; Nrange=3; NzBF=2; NBackFlowIdx=1;
  CHECK(BFDefIntCount()==44);
  storage[0]=storage[45]=123456789;
  BFBindDefTables(&cursor);
  CHECK(cursor==storage+45);
  CHECK(read_range(3,0,-1,NULL,NULL)==0);
  memcpy(expected,storage+1,sizeof(expected));
  CHECK(read_range(4,0,-1,NULL,NULL)==0);
  CHECK(memcmp(expected,storage+1,sizeof(expected))==0);
  CHECK(BFSeamPhase[0][2]==0);
  CHECK(read_range(3,1,-1,NULL,NULL)!=0);
  CHECK(read_range(4,0,2,"0 3 1 -1\n",NULL)!=0);
  CHECK(read_range(4,1,0,"0 0 0 -1\n",NULL)!=0);
  CHECK(read_range(4,1,2,"0 3 1 1\n",NULL)!=0);
  for (i=0;i<sizeof(bad)/sizeof(bad[0]);i++) CHECK(read_range(4,1,1,bad[i],NULL)!=0);
  CHECK(read_range(4,1,-1,NULL,"\n")!=0);
  CHECK(read_range(4,1,-1,NULL,"0 0 0 1\n")!=0);
  {
    const char embedded[]="0 0 0 1\0junk\n";
    char longline[512];
    int values[4], columns;
    FILE *fp=input(embedded,sizeof(embedded)-1);
    CHECK(BFReadRangeRow(fp,values,&columns)!=0); fclose(fp);
    memset(longline,' ',sizeof(longline)); memcpy(longline,"0 0 0 1",7);
    fp=input(longline,sizeof(longline));
    CHECK(BFReadRangeRow(fp,values,&columns)!=0); fclose(fp);
    fp=input("0 0 0 1",7);
    CHECK(BFReadRangeRow(fp,values,&columns)==0 && columns==4); fclose(fp);
  }
  CHECK(read_range(4,1,-1,NULL,NULL)==0);
  QPTrans=trans; QPTransSgn=signs; NMPTrans=1; NQPTrans=2;
  CHECK(BFValidateSeamTransforms()==0);
  s1[3]=1; /* unused row does not impose an extra symmetry */
  CHECK(BFValidateSeamTransforms()==0);
  NMPTrans=2;
  CHECK(BFValidateSeamTransforms()!=0);
  s1[3]=-1;
  CHECK(BFValidateSeamTransforms()==0);
  NMPTrans=1; QPTrans=trans+1; QPTransSgn=signs+1;
  CHECK(BFValidateSeamTransforms()==0); /* even a single nonidentity row is checked */
  t1[1]=3; t1[2]=2;
  CHECK(BFValidateSeamTransforms()!=0); /* permutation maps a bond outside range */
  t1[1]=2; t1[2]=3;
  QPTrans=trans; QPTransSgn=signs; NMPTrans=2;
  iFlgOrbitalGeneral=1; NQPOptTrans=1; QPOptTrans=opt; QPOptTransSgn=optsigns;
  CHECK(BFValidateSeamTransforms()==0);
  s1[3]=1;
  CHECK(BFValidateSeamTransforms()!=0);
  CHECK(storage[0]==123456789 && storage[45]==123456789);
  BFFreeDefTables();
  CHECK(PosBF==NULL && RangeIdx==NULL && BFSeamPhase==NULL && BackFlowIdx==NULL);
  CHECK(remove("tmp/backflow_seam_unit.input")==0);
  if (created) CHECK(rmdir("tmp")==0);
  puts("BackFlow seam input contract passed");
  return 0;
}
