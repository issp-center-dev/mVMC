#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "version.h"

void StdFace_main(char *fname);

int main(int argc, char *argv[])
{
  if (strcmp(argv[1], "-v") == 0) {
    printVersion();
    exit(0);
  }
  else {
    StdFace_main(argv[1]);
  }
}
