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
#ifndef MVMC_VERSION_H
#define MVMC_VERSION_H

/* Semantic Versioning http://semver.org */
/* <major>.<minor>.<patch>-<prerelease> */
/*
The version number is defined only here.
CMakeLists.txt, dist.sh and doc/(en|ja)/source/conf.py read the four lines
below, so keep the form "#define NAME value".
*/
#define MVMC_VERSION_MAJOR  1
#define MVMC_VERSION_MINOR  4
#define MVMC_VERSION_PATCH  0
#define MVMC_VERSION_PRERELEASE  "" /* "alpha", "beta.1", etc. */

const char *MVMC_GetGitHash(void);
void MVMC_PrintVersion(void);

#endif /* MVMC_VERSION_H */
