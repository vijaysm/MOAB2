
#include "MBAlloc.hpp"
#include <stdio.h>


// C malloc/free fucntions
void* mdb_malloc(size_t size, const char* filename=NULL, int linenumber=0)
{
  filename = filename;
  linenumber = linenumber;
  return malloc(size);
}

void* mdb_calloc(size_t size, const char* filename=NULL, int linenumber=0)
{
  filename = filename;
  linenumber = linenumber;
  void* ptr = malloc(size);
  memset(ptr, 0, size);
  return ptr;
}

void mdb_free(void* ptr)
{
  free(ptr);
}



