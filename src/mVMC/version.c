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
#include <stdio.h>
#include <string.h>
#include "version.h"

/* version_git.h is written into the build directory by cmake/git_hash.cmake. */
#ifdef MVMC_HAVE_VERSION_GIT_H
#include "version_git.h"
#endif
#ifndef MVMC_GIT_HASH
#define MVMC_GIT_HASH ""
#endif

/* Abbreviated hash (8 digits) of the commit which mVMC was built from,
   followed by "-dirty" if the source had changes which were not committed.
   An empty string if the commit is not known. */
const char *MVMC_GetGitHash(void) {
  return MVMC_GIT_HASH;
}

void MVMC_PrintVersion(void) {
  printf("mVMC version %d.%d.%d",MVMC_VERSION_MAJOR,MVMC_VERSION_MINOR,MVMC_VERSION_PATCH);
  if(strlen(MVMC_VERSION_PRERELEASE)>0) {
    printf("-%s",MVMC_VERSION_PRERELEASE);
  }
  if(strlen(MVMC_GetGitHash())>0) {
    printf(" (%s)",MVMC_GetGitHash());
  }
  printf("\n");
  return;
}
