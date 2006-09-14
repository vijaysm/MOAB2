/*
 * Library for determining the type of a geometric model file
 *
 * Copyright 2006, Jason Kraftcheck (kraftche@cae.wisc.edu)
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "geom_file_type.h"

const char* get_geom_file_type( const char* name )
{
  FILE* file;
  const char* result = 0;
  
  file = fopen( name, "r" );
  if (file) {
    result = get_geom_fptr_type( file );
    fclose( file );
  }
  
  return result;
}

const char* get_geom_fptr_type( FILE* file )
{
  static const char* CUBIT_NAME = GF_CUBIT_FILE_TYPE;
  static const char*  STEP_NAME = GF_STEP_FILE_TYPE;
  static const char*  IGES_NAME = GF_IGES_FILE_TYPE;
  static const char*   SAT_NAME = GF_ACIS_TXT_FILE_TYPE;
  static const char*   SAB_NAME = GF_ACIS_BIN_FILE_TYPE;
  
  if (is_cubit_file(file))
    return CUBIT_NAME;
  else if (is_step_file(file))
    return STEP_NAME;
  else if (is_iges_file(file))
    return IGES_NAME;
  else if (is_acis_bin_file(file))
    return SAB_NAME;
  else if (is_acis_txt_file(file))
    return SAT_NAME;
  else
    return 0;
}

int is_cubit_file( FILE* file )
{
  unsigned char buffer[4];
  return !fseek(file, 0, SEEK_SET) &&
         fread(buffer, 4, 1, file) &&
         !memcmp(buffer, "CUBE", 4);
}

int is_step_file( FILE* file )
{
  unsigned char buffer[9];
  return !fseek(file, 0, SEEK_SET) &&
         fread(buffer, 9, 1, file) &&
         !memcmp(buffer, "ISO-10303", 9);
}

int is_iges_file( FILE* file )
{
  unsigned char buffer[10];
  return !fseek(file, 72, SEEK_SET) &&
         fread(buffer, 10, 1, file) &&
         !memcmp(buffer, "S      1\r\n", 10);
}

int is_acis_bin_file( FILE* file )
{
  char buffer[15];
  return !fseek(file, 0, SEEK_SET) &&
         fread(buffer, 15, 1, file) &&
         !memcmp(buffer, "ACIS BinaryFile", 9);
}

int is_acis_txt_file( FILE* file )
{
  char buffer[5];
  int version, length;
  
  if (fseek(file,0,SEEK_SET) || 
      2 != fscanf( file, "%d %*d %*d %*d %d ", &version, &length ))
    return 0;
    
  if (version < 1 || version >0xFFFF)
    return 0;
  
    // Skip appliation name
  if (fseek(file, length, SEEK_CUR))
    return 0;
    
    // Read length of version string followed by first 5 characters
  if (2 != fscanf(file, "%d %4s", &length, buffer))
    return 0;
    
  return !strcmp( buffer, "ACIS" );
}

#ifdef DEBUG
int main( int argc, char* argv[] )
{
  while (--argc) 
    printf("%20s : %s\n", argv[argc], get_geom_file_type(argv[argc]));
  return 0;
}
#endif
