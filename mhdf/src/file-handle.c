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

#include <stdlib.h>
#include <string.h>
#include <H5Ipublic.h>
#include "file-handle.h"
#include "status.h"
#include "util.h"

#define FILE_HANDLE_MAGIC 0xFEEDFEED

int mhdf_check_valid_file( FileHandle* handle, mhdf_Status* status )
{
  if (!handle)
  {
    mhdf_setFail( status, "NULL file handle." );
    return 0;
  }
  
  if (handle->magic !=  FILE_HANDLE_MAGIC)
  {
    mhdf_setFail( status, "Invalid file handle." );
    return 0;
  }
  
  return 1;
}

FileHandle* mhdf_alloc_FileHandle( hid_t hdf_table, mhdf_Status* status )
{
  FileHandle* rval = (FileHandle*)mhdf_malloc( sizeof(FileHandle), status );
  if (!rval) return NULL;
  
  rval->magic = FILE_HANDLE_MAGIC;
  rval->hdf_handle = hdf_table;
  rval->open_handle_count = 0;
  rval->elem_type_handles = 0;
  rval->num_type_handles = 0;
  rval->type_array_len = 0;
  rval->max_id = 0L;
  return rval;
}

int mhdf_add_elem_type( FileHandle* handle, hid_t hdf_handle, mhdf_Status* status )
{
  hid_t* new_array;
  size_t new_size;
  int result = handle->num_type_handles;
  
  if (handle->type_array_len == handle->num_type_handles)
  {
    new_size = handle->type_array_len * 2;
    if (handle->type_array_len == 0)
      new_size = 8;
    new_array = (hid_t*)mhdf_malloc( new_size * sizeof(hid_t), status );
    if (!new_array) return 0;
    
    if (handle->num_type_handles)
    {
      memcpy( new_array, handle->elem_type_handles, handle->num_type_handles );
      free( handle->elem_type_handles );
    }
    handle->type_array_len = new_size;
    handle->elem_type_handles = new_array;
  }
  
  handle->num_type_handles++;
  handle->elem_type_handles[result] = hdf_handle;
  return result;
}

hid_t mhdf_handle_from_type_index( FileHandle* handle, 
                              int index,
                              mhdf_Status* status )
{
  if ((unsigned int)index < (unsigned int)handle->num_type_handles)
    return handle->elem_type_handles[index];

  mhdf_setFail( status, "Invalid element type index." );
  return (hid_t)-1;
}
