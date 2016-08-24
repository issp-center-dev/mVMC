/*-------------------------------------------------------------
 * Variational Monte Carlo
 * split loop
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/

inline void SplitLoop(int *start, int *end,
                      const int loopLength, 
                      const int mpiRank, const int mpiSize) {
  int ist,ien;
  int idiv,imod;

  if(mpiSize<loopLength){
    imod = loopLength%mpiSize;
    if(imod==0){
      idiv = loopLength/mpiSize;
      ist = idiv*mpiRank;
      ien = ist+idiv;
    } else {
      idiv = (loopLength-imod)/mpiSize;
      if(mpiRank<mpiSize-imod) {
        ist = idiv*mpiRank;
        ien = ist+idiv;
      } else {
        ist = idiv*mpiRank + mpiRank-(mpiSize-imod);
        ien = ist+idiv+1;
      }
    }
  } else {
    if(mpiRank<loopLength) {
      ist = mpiRank;
      ien = mpiRank+1;
    } else {
      ist = loopLength;
      ien = loopLength;
    }
  }

  /* printf("loop=%d,rank=%d,size=%d,ist=%d,ien=%d\n", */
  /*        loopLength,mpiRank,mpiSize,ist,ien); */

  *start = ist;
  *end = ien;

  return;
}
