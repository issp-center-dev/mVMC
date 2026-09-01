#ifndef _GC_SIZE_H
#define _GC_SIZE_H

#include <stddef.h>

size_t GCCheckedMulSize(size_t a, size_t b);
size_t GCCheckedAddSize(size_t a, size_t b);
size_t GCNonnegativeSize(long value, const char *name);
void *GCCheckedMalloc(size_t nbytes, const char *name);
int GCCheckedSizeToInt(size_t value, const char *name);

#endif
