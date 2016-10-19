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
