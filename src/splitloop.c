/*
mVMC - A numerical solver package for a wide range of quantum lattice models based on many-variable Variational Monte Carlo method
Copyright (C) 2016 Takahiro Misawa, Satoshi Morita, Takahiro Ohgoe, Kota Ido, Mitsuaki Kawamura, Takeo Kato, Masatoshi Imada.

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
 * Variational Monte Carlo
 * split loop
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/
#ifndef _SRC_SPLITLOOP
#define _SRC_SPLITLOOP
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
#endif