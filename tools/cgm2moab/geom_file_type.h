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

#ifndef GEOM_FILE_TYPE_H
#define GEOM_FILE_TYPE_H

#ifdef __cplusplus
extern "C" {
#endif

#define GF_CUBIT_FILE_TYPE    "CUBIT"
#define GF_STEP_FILE_TYPE     "STEP"
#define GF_IGES_FILE_TYPE     "IGES"
#define GF_ACIS_TXT_FILE_TYPE "ACIS_SAT"
#define GF_ACIS_BIN_FILE_TYPE "ACIS_SAB"

/* Get the type of a file.
   Return value is one of the above constants
 */
const char* get_geom_file_type( const char* filename );
const char* get_geom_fptr_type( FILE* file );

int is_cubit_file( FILE* file );
int is_step_file( FILE* file );
int is_iges_file( FILE* file );
int is_acis_txt_file( FILE* file );
int is_acis_bin_file( FILE* file );

#ifdef __cplusplus
} // extern "C"
#endif

#endif
