/*
mVMC - A numerical solver package for a wide range of quantum lattice models based on many-variable Variational Monte Carlo method
Copyright (C) 2016 Takahiro Misawa, Satoshi Morita, Takahiro Ohgoe, Kota Ido, Mitsuaki Kawamura, Takeo Kato, Masatoshi Imada.

his program is developed based on the mVMC-mini program
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
#ifndef _INCLUDE_VERSION
#define _INCLUDE_VERSION

#include <stdio.h>
#include <string.h>

/* Semantic Versioning http://semver.org */
/* <major>.<minor>.<patch>-<prerelease> */
#define VERSION_MAJOR  0
#define VERSION_MINOR  0
#define VERSION_PATCH  2
#define VERSION_PRERELEASE  "" /* "alpha", "beta.1", etc. */


void printVersion() {
  printf("mVMC version %d.%d.%d",VERSION_MAJOR,VERSION_MINOR,VERSION_PATCH);
  if(strlen(VERSION_PRERELEASE)>0) {
    printf("-%s",VERSION_PRERELEASE);
  }
  printf("\n");
  return;
};

#endif /* _INCLUDE_VERSIN */
