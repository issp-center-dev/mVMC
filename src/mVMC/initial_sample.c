/*
mVMC - A numerical solver package for a wide range of quantum lattice models based on many-variable Variational Monte Carlo method
Copyright (C) 2016 The University of Tokyo, All rights reserved.

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
 * local spins of the initial sample
 *-------------------------------------------------------------*/
#include "global.h"
#include "initial_sample.h"

/* Initialize eleIdx, eleCfg (and eleSpn for FSZ) and place the local spins.
   The itinerant electrons are placed by the caller afterwards.
   eleSpn==NULL : the electron mi+si*Ne has spin si.
   eleSpn!=NULL : FSZ. The first Ne+Sz electrons are up and the rest are down.

   The direction of each local spin is drawn at random. Then the itinerant
   electrons of one spin may outnumber the sites without local spin, and the
   caller would look for a free site forever. In that case the local spins
   are drawn again, up to nTryMax times, and after that with the number of up
   local spins restricted to the values with which the electrons fit.
   Random numbers are consumed in the same way as before as long as the
   first draw fits.

   Return 0, or 1 if the electrons do not fit with any direction of the
   local spins. */
int InitialSamplePlaceLocalSpinTry(int *eleIdx, int *eleCfg, int *eleSpn,
                                   const int nTryMax) {
  const int nsize = Nsize;
  const int nsite2 = Nsite2;
  int ri,mi,si,msi,rsi;
  int nLoc=0,nFree,nUpMin,nUpMax,nUp,nRest;
  int nTotal[2],nLocSpin[2],nTry;
  int fit=0;

  if(eleSpn==NULL) {
    nTotal[0] = Ne;
  } else {
    /* TwoSz==-1: Sz is not conserved but initially we take Sz=0 */
    nTotal[0] = Ne + ((TwoSz==-1) ? 0 : TwoSz/2);
  }
  nTotal[1] = nsize - nTotal[0];
  for(ri=0;ri<Nsite;ri++) {
    if(LocSpn[ri]==1) nLoc++;
  }
  nFree = Nsite - nLoc;

  /* range of the number of up local spins with which the electrons fit */
  nUpMin = 0;
  if(nUpMin < nLoc-nTotal[1]) nUpMin = nLoc-nTotal[1];
  if(nUpMin < nTotal[0]-nFree) nUpMin = nTotal[0]-nFree;
  nUpMax = nLoc;
  if(nUpMax > nTotal[0]) nUpMax = nTotal[0];
  if(nUpMax > nFree-nTotal[1]+nLoc) nUpMax = nFree-nTotal[1]+nLoc;
  if(nUpMin > nUpMax) return 1;

  for(nTry=0;nTry<=nTryMax;nTry++) {
    /* initialize */
    #pragma omp parallel for default(shared) private(msi)
    for(msi=0;msi<nsize;msi++) eleIdx[msi] = -1;
    #pragma omp parallel for default(shared) private(rsi)
    for(rsi=0;rsi<nsite2;rsi++) eleCfg[rsi] = -1;
    if(eleSpn!=NULL) {
      for(msi=0;msi<nsize;msi++) eleSpn[msi] = (msi<nTotal[0]) ? 0 : 1;
    }
    nLocSpin[0] = nLocSpin[1] = 0;

    if(nTry<nTryMax) {
      /* local spin */
      for(ri=0;ri<Nsite;ri++) {
        if(LocSpn[ri]==1) {
          if(eleSpn==NULL) {
            do {
              mi = gen_rand32()%Ne;
              si = (genrand_real2()<0.5) ? 0 : 1;
            } while(eleIdx[mi+si*Ne]!=-1);
            msi = mi+si*Ne;
          } else {
            do {
              msi = gen_rand32()%nsize;
            } while(eleIdx[msi]!=-1);
            si = eleSpn[msi];
          }
          eleCfg[ri+si*Nsite] = (eleSpn==NULL) ? mi : msi;
          eleIdx[msi] = ri;
          nLocSpin[si]++;
        }
      }
    } else {
      /* local spin, with the number of up spins in [nUpMin,nUpMax] */
      nUp = nUpMin + gen_rand32()%(nUpMax-nUpMin+1);
      nRest = nLoc;
      for(ri=0;ri<Nsite;ri++) {
        if(LocSpn[ri]==1) {
          si = (genrand_real2()*nRest<nUp) ? 0 : 1;
          if(si==0) nUp--;
          nRest--;
          do {
            mi = gen_rand32()%nTotal[si];
            msi = (si==0) ? mi : mi+nTotal[0];
          } while(eleIdx[msi]!=-1);
          eleCfg[ri+si*Nsite] = (eleSpn==NULL) ? mi : msi;
          eleIdx[msi] = ri;
          nLocSpin[si]++;
        }
      }
    }

    fit = (nTotal[0]-nLocSpin[0]<=nFree && nTotal[1]-nLocSpin[1]<=nFree);
    if(fit) break;
  }

  return 0;
}

int InitialSamplePlaceLocalSpin(int *eleIdx, int *eleCfg, int *eleSpn) {
  return InitialSamplePlaceLocalSpinTry(eleIdx,eleCfg,eleSpn,100);
}
