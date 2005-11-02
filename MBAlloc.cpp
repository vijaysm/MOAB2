/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */


#include "MBAlloc.hpp"
#include <stdio.h>


// C malloc/free fucntions
void* mdb_malloc(size_t size, const char* filename, int linenumber)
{
  filename = filename;
  linenumber = linenumber;
  return malloc(size);
}

void* mdb_calloc(size_t size, const char* filename, int linenumber)
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



