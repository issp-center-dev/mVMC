#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void StdFace_main(char *fname);

int main(int argc, char *argv[])
{
  if (strcmp(argv[1], "-v") == 0) {
    printf("\n Version 0.1\n\n");
    exit(0);
  }
  else {
    StdFace_main(argv[1]);
  }
}
